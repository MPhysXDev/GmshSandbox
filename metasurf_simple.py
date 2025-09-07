import gmsh, math, sys

# Simplified metasurface generator: build ONE elementary pattern, replicate by translation.
# Physical surface groups: black, white, disk, left, right, up, down, bottom. Volume: substrate.

# ---------------- Parameters -----------------
model_name = "meta_simple"
substrate_thickness = 0.01
period_x = 0.003
period_y = 0.003
num_pat_x = 3   # odd
num_pat_y = 3   # odd
circle_r = 0.0005
h_target = 1e-4

if num_pat_x % 2 == 0 or num_pat_y % 2 == 0:
    raise ValueError("num_pat_x / num_pat_y must be odd")
if circle_r >= min(period_x, period_y):
    raise ValueError("circle_r must be < period size")

pat_dx = 2*period_x
pat_dy = 2*period_y
ox = num_pat_x//2
oy = num_pat_y//2

# --------------- Helpers ---------------------
all_pts = set()

def P(x,y,z=0):
    tag = gmsh.model.occ.addPoint(x,y,z)
    all_pts.add(tag)
    return tag

def half_disk(cx, cy, r, side):
    d = gmsh.model.occ.addDisk(cx, cy, 0, r, r)
    L = 10*r
    if side == 'top':    rect = gmsh.model.occ.addRectangle(cx - L, cy - L, 0, 2*L, L)
    if side == 'bottom': rect = gmsh.model.occ.addRectangle(cx - L, cy,     0, 2*L, L)
    if side == 'left':   rect = gmsh.model.occ.addRectangle(cx,     cy - L, 0, L, 2*L)
    if side == 'right':  rect = gmsh.model.occ.addRectangle(cx - L, cy - L, 0, L, 2*L)
    cut = gmsh.model.occ.intersect([(2,d)], [(2,rect)], removeObject=True, removeTool=True)
    return cut[0][0][1]

def quarter_disk(cx, cy, r, quad):
    # quad strings: xpyp,xnyp,xpyn,xnyn (sign of (x-cx),(y-cy)) retained
    d = gmsh.model.occ.addDisk(cx, cy, 0, r, r)
    L = 10*r
    # build quadrant rectangle that keeps the inward part (toward center (0,0))
    # We'll use signs implicit in names
    if quad == 'xpyp':   rect = gmsh.model.occ.addRectangle(cx,     cy,     0, L, L)
    if quad == 'xnyp':   rect = gmsh.model.occ.addRectangle(cx-L,   cy,     0, L, L)
    if quad == 'xpyn':   rect = gmsh.model.occ.addRectangle(cx,     cy-L,   0, L, L)
    if quad == 'xnyn':   rect = gmsh.model.occ.addRectangle(cx-L,   cy-L,   0, L, L)
    cut = gmsh.model.occ.intersect([(2,d)], [(2,rect)], removeObject=True, removeTool=True)
    return cut[0][0][1]

# --------------- Build single elementary pattern (centered) ---------------
# Four rectangles (no boolean difference here; disks remain independent surfaces)
rects_local = {}  # (lx,ly)->tag
for lx in (0,1):
    for ly in (0,1):
        x0 = -period_x + lx*period_x
        y0 = -period_y + ly*period_y
        rects_local[(lx,ly)] = gmsh.model.occ.addRectangle(x0, y0, 0, period_x, period_y)

# Disks: 1 full, 4 half, 4 quarter
full_disk = gmsh.model.occ.addDisk(0,0,0,circle_r,circle_r)
half_disks = [
    half_disk(0,  period_y, circle_r, 'top'),
    half_disk(0, -period_y, circle_r, 'bottom'),
    half_disk(-period_x, 0, circle_r, 'left'),
    half_disk( period_x, 0, circle_r, 'right')
]
quarter_disks = [
    quarter_disk(-period_x, -period_y, circle_r, 'xpyp'),
    quarter_disk( period_x, -period_y, circle_r, 'xnyp'),
    quarter_disk(-period_x,  period_y, circle_r, 'xpyn'),
    quarter_disk( period_x,  period_y, circle_r, 'xnyn')
]
base_disk_surfs = [full_disk] + half_disks + quarter_disks

# --------------- Replicate by translation ---------------
all_rect_surfs = []  # (tag,gx,gy)
all_disk_surfs = []
for px in range(-ox, ox+1):
    for py in range(-oy, oy+1):
        tx = px * pat_dx
        ty = py * pat_dy
        # rectangles
        for (lx,ly), tag in rects_local.items():
            new = gmsh.model.occ.copy([(2, tag)])[0][1]
            gmsh.model.occ.translate([(2,new)], tx, ty, 0)
            gx = (px + ox)*2 + lx
            gy = (py + oy)*2 + ly
            all_rect_surfs.append((new, gx, gy))
        # disks (all share physical 'disk')
        for d in base_disk_surfs:
            nd = gmsh.model.occ.copy([(2,d)])[0][1]
            gmsh.model.occ.translate([(2,nd)], tx, ty, 0)
            all_disk_surfs.append(nd)

gmsh.model.occ.synchronize()

# --------------- Substrate box & fragment (imprint) ---------------
half_wx = num_pat_x * period_x
half_wy = num_pat_y * period_y
xmin, xmax = -half_wx, half_wx
ymin, ymax = -half_wy, half_wy
box = gmsh.model.occ.addBox(xmin, ymin, -substrate_thickness, 2*half_wx, 2*half_wy, substrate_thickness)

pattern_surfs_dimtags = [(2,t) for (t,_,_) in all_rect_surfs] + [(2,t) for t in all_disk_surfs]
# Fragment to carve top face
gmsh.model.occ.fragment([(3,box)], pattern_surfs_dimtags)
# DO NOT remove duplicates (would collapse partition)
gmsh.model.occ.synchronize()

# --------------- Classification after fragment ---------------
all_s2 = gmsh.model.getEntities(2)
# Collect z=0
z0_surfs = []
for _,s in all_s2:
    bb = gmsh.model.getBoundingBox(2,s)
    if abs(bb[2]) < 1e-12 and abs(bb[5]) < 1e-12:
        z0_surfs.append(s)

fullA = math.pi*circle_r**2
halfA = 0.5*fullA
quartA = 0.25*fullA
A_tol = 0.05*fullA

black = []
white = []
disk = []

for s in z0_surfs:
    A = gmsh.model.occ.getMass(2,s)
    if (abs(A-fullA) <= A_tol or abs(A-halfA) <= A_tol or abs(A-quartA) <= A_tol):
        disk.append(s)
        continue
    # rectangle piece (with holes already removed by fragment)
    bb = gmsh.model.getBoundingBox(2,s)
    cx = 0.5*(bb[0]+bb[3])
    cy = 0.5*(bb[1]+bb[4])
    gx = int(math.floor((cx - xmin)/period_x + 1e-12))
    gy = int(math.floor((cy - ymin)/period_y + 1e-12))
    if (gx + gy) % 2 == 0:
        black.append(s)
    else:
        white.append(s)

# --------------- Lateral / bottom surfaces ---------------
left = []; right = []; down = []; up = []; bottom = []
for _,s in all_s2:
    bb = gmsh.model.getBoundingBox(2,s)
    x0,y0,z0,x1,y1,z1 = bb
    if abs(z0 + substrate_thickness) < 1e-12 and abs(z1 + substrate_thickness) < 1e-12:
        bottom.append(s); continue
    if z0 < -1e-15 and z1 > 1e-15:  # vertical
        if abs(x0 - xmin) < 1e-12 and abs(x1 - xmin) < 1e-12: left.append(s)
        elif abs(x0 - xmax) < 1e-12 and abs(x1 - xmax) < 1e-12: right.append(s)
        elif abs(y0 - ymin) < 1e-12 and abs(y1 - ymin) < 1e-12: down.append(s)
        elif abs(y0 - ymax) < 1e-12 and abs(y1 - ymax) < 1e-12: up.append(s)

# --------------- Mesh sizing ---------------
for p in all_pts:
    gmsh.model.mesh.setSize([(0,p)], h_target)

# --------------- Physical groups ---------------
next_tag = 1

def PG(surfs,name):
    global next_tag
    if not surfs: return
    pg = gmsh.model.addPhysicalGroup(2, surfs, tag=next_tag)
    gmsh.model.setPhysicalName(2, pg, name)
    next_tag += 1

PG(black,'black')
PG(white,'white')
PG(disk,'disk')
PG(left,'left')
PG(right,'right')
PG(up,'up')
PG(down,'down')
PG(bottom,'bottom')

sub_pg = gmsh.model.addPhysicalGroup(3, [box], tag=100)
gmsh.model.setPhysicalName(3, sub_pg, 'substrate')

# --------------- Mesh & write ---------------
gmsh.model.mesh.generate(3)
gmsh.write('metasurface_supercell_simple.msh')

# Print quick counts
print(f"black={len(black)} white={len(white)} disk={len(disk)} lateral={len(left)+len(right)+len(up)+len(down)} bottom={len(bottom)}")

gmsh.finalize()
