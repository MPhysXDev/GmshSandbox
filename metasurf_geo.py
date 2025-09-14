# Author: Thierry Valet
# Copyright (C) 2025 - All rights reserved

"""
This script creates a chequerboard metasurface unit cell,
and then a supercell by replicating the unit cell in x and y directions.
The metasurface supercell is defined at the top surface of a substrate
of given thickness, with a source plane at a certain distance above the metasurface.
We use the native Gmsh kernel for the geometry definition.
A 3D mesh is generated for full-wave electromagnetic simulations.
"""

# Import necessary modules
import sys
import numpy as np
import gmsh

# Function to create a sub-unit cell with one rectangle and 
# arc-circle cutouts at each corner
def sub_unit_cell(
        x0, # Lower left corner x-coordinate
        y0, # Lower left corner y-coordinate
        wx, # Width of the rectangles
        wy, # Height of the rectangles
        r,  # Radius of the circular regions
        hx, # Mesh size at rectangle x edges
        hy, # Mesh size at rectangle y edges
        hc   # Mesh size at circular perimeter
    ):
    # Creates the lower left sector
    p0 = gmsh.model.geo.addPoint(x0, y0, 0, hc)
    p1 = gmsh.model.geo.addPoint(x0+r, y0, 0, hc)
    p2 = gmsh.model.geo.addPoint(x0, y0+r, 0, hc)

    l0 = gmsh.model.geo.addLine(p0, p1)
    c0 = gmsh.model.geo.addCircleArc(p1, p0, p2)
    l1 = gmsh.model.geo.addLine(p2, p0)
    ll0 = gmsh.model.geo.addCurveLoop([l0, c0, l1])
    s0 = gmsh.model.geo.addPlaneSurface([ll0])

    p3 = gmsh.model.geo.addPoint(x0+wx/2, y0, 0, hx)

    # Creates the lower right sector
    p4 = gmsh.model.geo.addPoint(x0+wx, y0, 0, hc)
    p5 = gmsh.model.geo.addPoint(x0+wx-r, y0, 0, hc)
    p6 = gmsh.model.geo.addPoint(x0+wx, y0+r, 0, hc)  

    l2 = gmsh.model.geo.addLine(p5, p4)
    c1 = gmsh.model.geo.addCircleArc(p6, p4, p5)
    l3 = gmsh.model.geo.addLine(p4, p6)
    ll1 = gmsh.model.geo.addCurveLoop([l2, l3, c1])
    s1 = gmsh.model.geo.addPlaneSurface([ll1])

    p7 = gmsh.model.geo.addPoint(x0+wx, y0+wy/2, 0, hy)

    # Creates the upper right sector
    p8 = gmsh.model.geo.addPoint(x0+wx, y0+wy, 0, hc)
    p9 = gmsh.model.geo.addPoint(x0+wx, y0+wy-r, 0, hc)
    p10 = gmsh.model.geo.addPoint(x0+wx-r, y0+wy, 0, hc)

    l4 = gmsh.model.geo.addLine(p9, p8)
    c2 = gmsh.model.geo.addCircleArc(p10, p8, p9)
    l5 = gmsh.model.geo.addLine(p8, p10)
    ll2 = gmsh.model.geo.addCurveLoop([l5, c2, l4])
    s2 = gmsh.model.geo.addPlaneSurface([ll2])

    p11 = gmsh.model.geo.addPoint(x0+wx/2, y0+wy, 0, hx)
    
    # Creates the upper left sector
    p12 = gmsh.model.geo.addPoint(x0, y0+wy, 0, hc)
    p13 = gmsh.model.geo.addPoint(x0, y0+wy-r, 0, hc)
    p14 = gmsh.model.geo.addPoint(x0+r, y0+wy, 0, hc)

    l6 = gmsh.model.geo.addLine(p13, p12)
    c3 = gmsh.model.geo.addCircleArc(p14, p12, p13)
    l7 = gmsh.model.geo.addLine(p12, p14)
    ll3 = gmsh.model.geo.addCurveLoop([-l6, -l7, -c3])
    s3 = gmsh.model.geo.addPlaneSurface([ll3])   

    p15 = gmsh.model.geo.addPoint(x0, y0+wy/2, 0, hy)

    l8 = gmsh.model.geo.addLine(p1, p3)
    l9 = gmsh.model.geo.addLine(p3, p5)
    l10 = gmsh.model.geo.addLine(p6, p7)
    l11 = gmsh.model.geo.addLine(p7, p9)
    l12 = gmsh.model.geo.addLine(p10, p11)
    l13 = gmsh.model.geo.addLine(p11, p14)
    l14 = gmsh.model.geo.addLine(p13, p15)
    l15 = gmsh.model.geo.addLine(p15, p2)
    ll4 = gmsh.model.geo.addCurveLoop([l8, l9, -c1, l10, l11,-c2,
                                       l12, l13, c3, l14, l15, -c0])
    s4 = gmsh.model.geo.addPlaneSurface([ll4])
    
    # List of created surfaces: lower left, lower right, 
    # upper right, upper left, center (rectangle with cutouts)
    surfaces = [s0, s1, s2, s3, s4]

    # List of edges
    front = [l0, l8, l9, l2]
    right = [l3, l10, l11, l4]
    back = [l6, -l13, -l12, -l5]
    left = [-l1, -l15, -l14, -l7]
    edges = [front, right, back, left]

    # List of corner points
    corners = [p0, p4, p8, p12]

    return [surfaces , edges, corners]

# Function to create a unit cell with four rectangle and 
# arc-circle cutouts at each inner and outer vertices
def unit_cell(
        x0, # Bottom left corner x-coordinate
        y0, # Bottom left corner y-coordinate
        wx, # Width of the rectangles
        wy, # Height of the rectangles
        r,  # Radius of the circular regions
        hx, # Mesh size at rectangle x edges
        hy, # Mesh size at rectangle y edges
        hc   # Mesh size at circular perimeter
    ):
    
    lower_left = sub_unit_cell(x0, y0, wx, wy, r, hx, hy, hc)
    lower_right = sub_unit_cell(x0+wx, y0, wx, wy, r, hx, hy, hc)
    upper_right = sub_unit_cell(x0+wx, y0+wy, wx, wy, r, hx, hy, hc)
    upper_left = sub_unit_cell(x0, y0+wy, wx, wy, r, hx, hy, hc)

    # arrays of surfaces in the unit cell
    black_rectangles = [lower_left[0][4], upper_right[0][4]]
    white_rectangles = [lower_right[0][4], upper_left[0][4]]
    cb = [
        [lower_left[0][0]], 
        [lower_left[0][1], lower_right[0][0]],
        [lower_right[0][1]]
    ]
    cm = [
        [lower_left[0][3], upper_left[0][0]],
        [lower_left[0][2], lower_right[0][3], upper_right[0][0], upper_left[0][1]],
        [lower_right[0][2], upper_right[0][1]]
    ]
    ct = [
        [upper_left[0][3]],
        [upper_left[0][2], upper_right[0][3]],
        [upper_right[0][2]]
    ]
    disks = [cb, cm, ct]

    # list of edges in the unit cell
    front = lower_left[1][0] + lower_right[1][0]
    right = lower_right[1][1] + upper_right[1][1]
    back = upper_left[1][2] + upper_right[1][2]
    left = lower_left[1][3] + upper_left[1][3]
    edges = [front, right, back, left]

    # list of corner points
    corners = [lower_left[2][0], lower_right[2][1], upper_right[2][2], upper_left[2][3]]

    return [black_rectangles, white_rectangles, disks, edges, corners]

# Function to create a full period of the metasurface by replicating the unit cell
def meta_period(
        wx, # Width of the rectangles
        wy, # Height of the rectangles
        r,  # Radius of the circular regions
        hx, # Mesh size at rectangle x edges
        hy, # Mesh size at rectangle y edges
        hc, # Mesh size at circular perimeter
        nx, # replication count in x
        ny  # replication count in y
    ):
    # Lists to store the different parts of the metasurface
    black = []
    white = []
    corner_disks = []
    x0_edge = []
    x1_edge = []
    y0_edge = []
    y1_edge = []
    inner_disks = []
    # Aggregate outer perimeter edges as chains of curves
    front_edges = []
    right_edges = []
    back_edges = []
    left_edges = []
    # List of corner point tags
    corner_points = []
    # Global node accumulator for quarter-circle surfaces per lattice node
    # Grid of size (2*ny + 1) x (2*nx + 1): each tile contributes to a 3x3 block
    inner_nodes = [[[] for _ in range(2*nx + 1)] for _ in range(2*ny + 1)]
    # Iterate over the number of replications in y and x 
    for j in range(ny):
        for i in range(nx):
            # Create the unit cell at the appropriate position
            out = unit_cell((2*i - nx) * wx, (2*j - ny) * wy, wx, wy, r, hx, hy, hc)
            # Assemble lists corresponding to the different parts of the metasurface
            # **********************************************************************
            # chequerboard rectangles
            black += out[0]
            white += out[1]
            # corner cutouts
            if ((i == 0) and (j == 0)):
                corner_disks += out[2][0][0]
                corner_points.append(out[4][0])
            if ((i == nx-1) and (j == 0)):
                corner_disks += out[2][0][2]
                corner_points.append(out[4][1])
            if ((i == 0) and (j == ny-1)):
                corner_disks += out[2][2][0]
                corner_points.append(out[4][3])
            if ((i == nx-1) and (j == ny-1)):
                corner_disks += out[2][2][2]
                corner_points.append(out[4][2])
            # x edge cutouts
            if (j == 0):
                if (i != 0):
                    x0_edge[-1] += out[2][0][0]
                x0_edge.append(out[2][0][1])
                if (i != nx-1):
                    x0_edge.append(out[2][0][2])
            if (j == ny-1):
                if (i != 0):
                    x1_edge[-1] += out[2][2][0]
                x1_edge.append(out[2][2][1])
                if (i != nx-1):
                    x1_edge.append(out[2][2][2])
            # y edge cutouts
            if (i == 0):
                if (j != 0):
                    y0_edge[-1] += out[2][0][0]
                y0_edge.append(out[2][1][0])
                if (j != ny-1):
                    y0_edge.append(out[2][2][0])
            if (i == nx-1):
                if (j != 0):
                    y1_edge[-1] += out[2][0][2]
                y1_edge.append(out[2][1][2])
                if (j != ny-1):
                    y1_edge.append(out[2][2][2])
            # collect outer perimeter edges
            if (j == 0):
                front_edges += out[3][0]
            if (i == nx-1):
                right_edges += out[3][1]
            if (j == ny-1):
                back_edges += out[3][2]
            if (i == 0):
                left_edges += out[3][3]
            # Accumulate all 3x3 cutout clusters from this tile into the global node grid
            # out[2] is a 3x3 array of lists: rows cb, cm, ct and columns left, center, right
            for rr in range(3):
                for cc in range(3):
                    cells = out[2][rr][cc]
                    if cells:
                        inner_nodes[2*j + rr][2*i + cc] += cells
    x_edge_disks =[a + b for a, b in zip(x0_edge, x1_edge)]
    y_edge_disks =[a + b for a, b in zip(y0_edge, y1_edge)]

    # Collect all interior node cuts (exclude outer boundary nodes)
    # Order: left-to-right within each row, bottom-to-top across rows
    for jn in range(1, 2*ny):
        for in_idx in range(1, 2*nx):
            cuts = inner_nodes[jn][in_idx]
            if cuts:
                inner_disks.append(cuts)

    return ( 
        black, white, corner_disks, x_edge_disks, y_edge_disks, inner_disks, 
        front_edges, right_edges, back_edges, left_edges,
        corner_points
    )

"""
def meta_model(wx, wy, r, hx, hy, hc, hb, nx, ny, thickness):
    # Build metasurface period and gather outer chains
    (
        black, white, corners, x_edge, y_edge, inner_array,
        front_edges, right_edges, back_edges, left_edges,
    ) = meta_period(wx, wy, r, hx, hy, hc, nx, ny)

    # With known orientations: front (L->R), right (F->B), back (L->R), left (F->B)
    # Get the top corner point tags from chain endpoints
    fl_top, fr_top = chain_start_end(front_edges)
    fr_top2, br_top = chain_start_end(right_edges)
    bl_top, br_top2 = chain_start_end(back_edges)
    fl_top2, bl_top2 = chain_start_end(left_edges)

    # Use chains directly (no reorientation)
    front_chain = front_edges
    right_chain = right_edges
    back_chain = back_edges
    left_chain = left_edges

    # Domain extents for bottom rectangle
    x_min, x_max = -nx*wx, nx*wx
    y_min, y_max = -ny*wy, ny*wy

    # Create bottom rectangle corner points (z = -thickness)
    p_bl = gmsh.model.geo.addPoint(x_min, y_min, -thickness, hb)
    p_br = gmsh.model.geo.addPoint(x_max, y_min, -thickness, hb)
    p_tr = gmsh.model.geo.addPoint(x_max, y_max, -thickness, hb)
    p_tl = gmsh.model.geo.addPoint(x_min, y_max, -thickness, hb)

    # Bottom perimeter straight lines
    e_bottom = gmsh.model.geo.addLine(p_bl, p_br)
    e_right = gmsh.model.geo.addLine(p_br, p_tr)
    e_top = gmsh.model.geo.addLine(p_tr, p_tl)
    e_left = gmsh.model.geo.addLine(p_tl, p_bl)

    # Create vertical corner edges from top endpoints to bottom corners
    # If some endpoints are None (degenerate), skip side surface generation
    v_fl = add_vertical(fl_top, p_bl)
    v_fr = add_vertical(fr_top, p_br)
    v_br = add_vertical(br_top, p_tr)
    v_bl = add_vertical(bl_top, p_tl)

    # Build side surfaces using curve loops: top chain + vertical + bottom edge (reversed) + vertical (reversed)
    side_front = None
    side_right = None
    side_back = None
    side_left = None
    if all(v is not None for v in [v_fl, v_fr]):
        loop_front = gmsh.model.geo.addCurveLoop(front_chain + [v_fr, -e_bottom, -v_fl])
        side_front = gmsh.model.geo.addSurfaceFilling([loop_front])
    if all(v is not None for v in [v_fr, v_br]):
        loop_right = gmsh.model.geo.addCurveLoop(right_chain + [v_br, -e_right, -v_fr])
        side_right = gmsh.model.geo.addSurfaceFilling([loop_right])
    if all(v is not None for v in [v_bl, v_br]):
        loop_back = gmsh.model.geo.addCurveLoop(back_chain + [v_br, e_top, -v_bl])
        side_back = gmsh.model.geo.addSurfaceFilling([loop_back])
    if all(v is not None for v in [v_fl, v_bl]):
        loop_left = gmsh.model.geo.addCurveLoop(left_chain + [v_bl, e_left, -v_fl])
        side_left = gmsh.model.geo.addSurfaceFilling([loop_left])

    # Bottom surface
    loop_bottom = gmsh.model.geo.addCurveLoop([e_bottom, e_right, e_top, e_left])
    bottom_surface = gmsh.model.geo.addPlaneSurface([loop_bottom])

    # Collect all top surfaces to close the volume
    top_surfaces = []
    top_surfaces += black
    top_surfaces += white
    top_surfaces += corners
    for arr in x_edge:
        top_surfaces += arr
    for arr in y_edge:
        top_surfaces += arr
    for arr in inner_array:
        top_surfaces += arr

    # Create volume from surface loop
    shell_surfaces = [bottom_surface]
    if side_front is not None:
        shell_surfaces.append(side_front)
    if side_right is not None:
        shell_surfaces.append(side_right)
    if side_back is not None:
        shell_surfaces.append(side_back)
    if side_left is not None:
        shell_surfaces.append(side_left)
    shell_surfaces += top_surfaces

    sl = gmsh.model.geo.addSurfaceLoop(shell_surfaces)
    vol = gmsh.model.geo.addVolume([sl])

    substrate = {
        'volume': vol,
        'bottom': bottom_surface,  # bottom plane at z = -thickness
        'front': [side_front] if side_front is not None else [],  # y-min side
        'back': [side_back] if side_back is not None else [],     # y-max side
        'left': [side_left] if side_left is not None else [],
        'right': [side_right] if side_right is not None else [],
    }

    return black, white, corners, x_edge, y_edge, inner_array, substrate
"""

# Configuration parameters
model_name = 'meta_surf'
num_x = 5          # replication count in x
num_y = 5          # replication count in y
rect_x = 0.003     # rectangle dimension along x
rect_y = 0.003     # rectangle dimension along y
circle_r = 0.0005  # disk radius
hx = 1e-3          # mesh size at rectangle x edges
hy = 1e-3          # mesh size at rectangle y edges
hc = 5e-4          # mesh size at circular perimeter
hb = 5e-3          # mesh size at bottom rectangle corner points
substrate_thickness = 0.01
source_distance = 0.06
absorber_distance = 0.06

# Initialize Gmsh and create the model
gmsh.initialize()
gmsh.model.add(model_name)

# Generate period and build substrate volume via geometric construction (no extrusion)
( black, white, corner_disks, x_edge_disks, y_edge_disks, inner_disks, 
        front_edges, right_edges, back_edges, left_edges,
        corner_points ) = meta_period(rect_x, rect_y, circle_r, hx, hy, hc, num_x, num_y)

# Remove duplicate entities
gmsh.model.geo.removeAllDuplicates()

# Synchronize the CAD kernel with the Gmsh model
gmsh.model.geo.synchronize()

# Define physical groups for the metasurface components
gmsh.model.addPhysicalGroup(2, black, tag=1, name="black_rectangles")
gmsh.model.addPhysicalGroup(2, white, tag=2, name="white_rectangles")
gmsh.model.addPhysicalGroup(2, corner_disks, tag=10, name="corner_disks")
for i, cut in enumerate(x_edge_disks): 
    gmsh.model.addPhysicalGroup(2, cut, tag=100+i, name=f"x_edge_disk_{i}")
for i, cut in enumerate(y_edge_disks): 
    gmsh.model.addPhysicalGroup(2, cut, tag=200+i, name=f"y_edge_disk_{i}")
for i, cut in enumerate(inner_disks):
    gmsh.model.addPhysicalGroup(2, cut, tag=1000+i, name=f"inner_disk_{i}")

"""
# New physical groups for the substrate (renamed faces)
if substrate.get('volume') is not None:
    gmsh.model.addPhysicalGroup(3, [substrate['volume']], tag=5000, name="substrate")
if substrate.get('bottom') is not None:
    gmsh.model.addPhysicalGroup(2, [substrate['bottom']], tag=10000, name="substrate_bottom")
if substrate.get('front'):
    gmsh.model.addPhysicalGroup(2, substrate['front'], tag=20000, name="substrate_front")
if substrate.get('back'):
    gmsh.model.addPhysicalGroup(2, substrate['back'], tag=30000, name="substrate_back")
if substrate.get('left'):
    gmsh.model.addPhysicalGroup(2, substrate['left'], tag=40000, name="substrate_left")
if substrate.get('right'):
    gmsh.model.addPhysicalGroup(2, substrate['right'], tag=50000, name="substrate_right")
"""

# Display the generated metasurface
gmsh.fltk.run()

gmsh.finalize()