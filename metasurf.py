import gmsh
import math
import sys

"""
Metasurface supercell generator (Option 1: literal replication of elementary pattern)

Elementary pattern: 2x2 rectangles (each period_x x period_y) at z=0 forming
footprint 2*period_x by 2*period_y, centered at pattern center. Pattern includes:
 - 1 full disk at internal junction (pattern center)
 - 4 half disks at midpoints of outer edges
 - 4 quarter disks at outer corners
All disks (full/half/quarter) share the physical tag 'disk'. Rectangles alternate
black/white via (global_rect_i + global_rect_j) % 2.

Patterns are replicated on lattice (2*period_x, 2*period_y) for an odd number of
patterns num_cells_x, num_cells_y, centered so central pattern center at (0,0).

Substrate: single box volume spanning full supercell in xy, thickness substrate_thickness
extending from z=-substrate_thickness to z=0. Pattern surfaces imprint (fragment)
the top surface of substrate. Lateral faces tagged left/right/down/up, bottom tagged bottom.

Mesh: target_mesh_size set on all key points (rectangle corners, disk centers, arc endpoints).
Periodic-readiness: we also enforce consistent subdivisions on outer boundary curves via
setTransfiniteCurve so opposite sides have identical discretization counts.
"""

config = {
    'model_name': 'meta0',
    'geometry': {
        'substrate_thickness': 0.01,
        'period_x': 0.003,
        'period_y': 0.003,
        'num_cells_x': 3,      # must be odd
        'num_cells_y': 3,      # must be odd
        'circle_radius': 0.0005,
        'target_mesh_size': 1e-4
    }
}

g = config['geometry']
"""Simplified metasurface generator (minimal) with optional interactive visualization.
Flags (simple parsing):
  --preview : open GUI after geometry (before meshing)
  --gui     : open GUI after meshing
"""

import gmsh, math, sys

# Configuration block (single source of parameters)
config = {
    'model_name': 'meta_simplified',
    'substrate_thickness': 0.01,
    'period_x': 0.003,
    'period_y': 0.003,
    'num_pat_x': 3,      # must be odd
    'num_pat_y': 3,      # must be odd
    'circle_r': 0.0005,
    'mesh_size': 1e-4
}

model_name = config['model_name']
substrate_thickness = config['substrate_thickness']
period_x = config['period_x']
period_y = config['period_y']
num_pat_x = config['num_pat_x']
num_pat_y = config['num_pat_y']
circle_r = config['circle_r']
h_target = config['mesh_size']

# ---------------- Parameters (edit here) -----------------

# No CLI flags; always show GUI before finalize.

if num_pat_x % 2 == 0 or num_pat_y % 2 == 0:
    raise ValueError("num_pat_x / num_pat_y must be odd")
if circle_r >= min(period_x, period_y):
    raise ValueError("circle_r must be < period size")

gmsh.initialize()
gmsh.model.add(model_name)

pat_dx = 2*period_x
pat_dy = 2*period_y
ox = num_pat_x//2
oy = num_pat_y//2

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
    d = gmsh.model.occ.addDisk(cx, cy, 0, r, r)
    L = 10*r
    if quad == 'xpyp':   rect = gmsh.model.occ.addRectangle(cx,   cy,   0, L, L)
    if quad == 'xnyp':   rect = gmsh.model.occ.addRectangle(cx-L, cy,   0, L, L)
    if quad == 'xpyn':   rect = gmsh.model.occ.addRectangle(cx,   cy-L, 0, L, L)
    if quad == 'xnyn':   rect = gmsh.model.occ.addRectangle(cx-L, cy-L, 0, L, L)
    cut = gmsh.model.occ.intersect([(2,d)], [(2,rect)], removeObject=True, removeTool=True)
    return cut[0][0][1]

############################
# Elementary pattern build #
############################
base_rects = {}
for lx in (0,1):
    for ly in (0,1):
        x0 = -period_x + lx*period_x
        y0 = -period_y + ly*period_y
        base_rects[(lx,ly)] = gmsh.model.occ.addRectangle(x0, y0, 0, period_x, period_y)

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
disk_base = [full_disk] + half_disks + quarter_disks

# Cut disks out of each rectangle ONCE so holes exist (retain disks)
cut_rects = {}
tools = [(2,d) for d in disk_base]
for key, rtag in base_rects.items():
    cut = gmsh.model.occ.cut([(2,rtag)], tools, removeObject=True, removeTool=False)
    cut_rects[key] = cut[0][0][1] if cut and cut[0] else rtag

###########
# Copying #
###########
all_rects = []  # (surface_tag, local_lx, local_ly)
all_disks = []
black = []
white = []
disk = []
for px in range(-ox, ox+1):
    for py in range(-oy, oy+1):
        tx = px * pat_dx; ty = py * pat_dy
        for (lx,ly), rtag in cut_rects.items():
            nt = gmsh.model.occ.copy([(2,rtag)])[0][1]
            gmsh.model.occ.translate([(2,nt)], tx, ty, 0)
            all_rects.append((nt,lx,ly))
            if (lx + ly) % 2 == 0:
                black.append(nt)
            else:
                white.append(nt)
        for d in disk_base:
            nd = gmsh.model.occ.copy([(2,d)])[0][1]
            gmsh.model.occ.translate([(2,nd)], tx, ty, 0)
            all_disks.append(nd)
            disk.append(nd)

gmsh.model.occ.synchronize()

# Substrate & fragment
half_wx = num_pat_x * period_x
half_wy = num_pat_y * period_y
xmin,xmax = -half_wx, half_wx
ymin,ymax = -half_wy, half_wy
box = gmsh.model.occ.addBox(xmin, ymin, -substrate_thickness, 2*half_wx, 2*half_wy, substrate_thickness)
pat_dimtags = [(2,t) for (t,_,_) in all_rects] + [(2,t) for t in all_disks]
gmsh.model.occ.fragment([(3,box)], pat_dimtags)
gmsh.model.occ.synchronize()

# (Removed preview phase) Build then mesh, then display once at end.

surf2 = gmsh.model.getEntities(2)

left=[]; right=[]; up=[]; down=[]; bottom=[]
for _,s in surf2:
    bb = gmsh.model.getBoundingBox(2,s)
    x0,y0,z0,x1,y1,z1 = bb
    if abs(z0 + substrate_thickness)<1e-12 and abs(z1 + substrate_thickness)<1e-12:
        bottom.append(s); continue
    if z0 < -1e-15 and z1 > 1e-15:
        if abs(x0 - xmin)<1e-12 and abs(x1 - xmin)<1e-12: left.append(s)
        elif abs(x0 - xmax)<1e-12 and abs(x1 - xmax)<1e-12: right.append(s)
        elif abs(y0 - ymin)<1e-12 and abs(y1 - ymin)<1e-12: down.append(s)
        elif abs(y0 - ymax)<1e-12 and abs(y1 - ymax)<1e-12: up.append(s)

# Physical groups with fixed tags
def add_pg(tag, surfs, name):
    if not surfs: return
    gmsh.model.addPhysicalGroup(2, surfs, tag)
    gmsh.model.setPhysicalName(2, tag, name)

add_pg(1, black, 'black')
add_pg(2, white, 'white')
add_pg(3, disk,  'disk')
add_pg(4, left,  'left')
add_pg(5, right, 'right')
add_pg(6, up,    'up')
add_pg(7, down,  'down')
add_pg(8, bottom,'bottom')
gmsh.model.addPhysicalGroup(3, [box], 100)
gmsh.model.setPhysicalName(3, 100, 'substrate')

gmsh.model.mesh.generate(3)
gmsh.write(f"{model_name}.msh")
print(f"black={len(black)} white={len(white)} disk={len(disk)} left={len(left)} right={len(right)} up={len(up)} down={len(down)} bottom={len(bottom)}")

try:
    gmsh.fltk.run()
except Exception:
    pass

gmsh.finalize()
sys.exit(0)



