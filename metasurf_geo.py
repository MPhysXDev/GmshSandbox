# Author: Thierry Valet
# Copyright (C) 2025 - All rights reserved

"""
This script creates a chequerboard metasurface unit cell,
and then a supercell by replicating the unit cell in x and y directions.
The metasurface supercell is defined at the top surface of a substrate
of given thickness, with a source plane at a certain distance above the metasurface,
and an absorber plane at a further distance above the source plane.
We use the native Gmsh kernel for the geometry definition.
A 3D mesh is generated for full-wave electromagnetic simulations.
--- ALL DIMENSIONS IN METERS --- METASURFACE AT z=0, CENTERED AT (x,y) ORIGIN ---
--- SUBSTRATE BELOW, SOURCE AND ABSORBER ABOVE ---
"""

# Import necessary modules
import sys
import os
import gmsh

# General parameters
path_to_mesh = './'  # path to location to save the mesh file
display_model = True  # Set to True to display the model in Gmsh GUI
save_mesh = True  # Set to True to save the mesh to file

# Model parameters
model_name = 'meta_surf'
num_x = 10          # replication count in x
num_y = 10          # replication count in y
rect_x = 0.003     # rectangle dimension along x
rect_y = 0.003     # rectangle dimension along y
circle_r = 0.0005  # disk cutout radius
hx = 1e-3          # mesh size at rectangle x edges
hy = 1e-3          # mesh size at rectangle y edges
hc = 5e-4          # mesh size at circular cutout perimeters
hb = 5e-3          # mesh size at substrate bottom corners
hs = 5e-3          # mesh size at source plane corners
ha = 5e-3          # mesh size at absorber plane corners
substrate_thickness = 0.01
source_distance = 0.06      # distance from metasurface to source plane
absorber_distance = 0.12    # distance from metasurface to absorber plane

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
    back = [l7, -l13, -l12, -l5]
    left = [-l1, -l15, -l14, l6]
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

# Funtion to create the full 3D model
def meta_model(
        wx, # Width of the rectangles
        wy, # Height of the rectangles
        r,  # Radius of the circular regions
        nx, # Number of unit cells in x direction
        ny, # Number of unit cells in y direction
        subst_t, # Substrate thickness
        source_distance, # Distance from metasurface to source plane
        absorber_distance, # Distance from metasurface to absorber plane
        hx, # Mesh size at rectangle x edges
        hy, # Mesh size at rectangle y edges
        hc, # Mesh size at circular perimeter
        hb, # Mesh size at bottom rectangle corner points
        hs, # Mesh size at source plane corners
        ha # Mesh size at absorber plane corners
    ):

    # Create one period of the metasurface
    # ************************************
    (
        black, white, corner_disks, x_edge_disks, y_edge_disks, inner_disks,
        front_edges, right_edges, back_edges, left_edges,
        corner_points
    ) = meta_period(wx, wy, r, hx, hy, hc, nx, ny)

    # Remove duplicate entities
    gmsh.model.geo.removeAllDuplicates()

    # Create the substrate volume and faces
    # *************************************
    # Create bottom face
    x_min, x_max = -nx*wx, nx*wx 
    y_min, y_max = -ny*wy, ny*wy
    p_bfl = gmsh.model.geo.addPoint(x_min, y_min, -subst_t, hb)
    p_bfr = gmsh.model.geo.addPoint(x_max, y_min, -subst_t, hb)
    p_bbr = gmsh.model.geo.addPoint(x_max, y_max, -subst_t, hb)
    p_bbl = gmsh.model.geo.addPoint(x_min, y_max, -subst_t, hb)

    l_bbf = gmsh.model.geo.addLine(p_bfr, p_bfl)
    l_bbl = gmsh.model.geo.addLine(p_bfl, p_bbl)
    l_bbb = gmsh.model.geo.addLine(p_bbl, p_bbr)
    l_bbr = gmsh.model.geo.addLine(p_bbr, p_bfr)

    ll_bot = gmsh.model.geo.addCurveLoop([l_bbf, l_bbl, l_bbb, l_bbr])
    s_bot = gmsh.model.geo.addPlaneSurface([ll_bot])

    # Create front face
    l_fl = gmsh.model.geo.addLine(p_bfl, corner_points[0])
    l_fr = gmsh.model.geo.addLine(corner_points[1], p_bfr)

    ll_front = [l_bbf, l_fl] + front_edges + [l_fr]
    s_front = gmsh.model.geo.addPlaneSurface([gmsh.model.geo.addCurveLoop(ll_front)])

    # Create back face
    l_bl = gmsh.model.geo.addLine(p_bbl, corner_points[2])
    l_br = gmsh.model.geo.addLine(corner_points[3], p_bbr)

    ll_back = [-l_bbb, l_bl] + back_edges + [l_br]
    s_back = gmsh.model.geo.addPlaneSurface([gmsh.model.geo.addCurveLoop(ll_back)])

    # Create left face
    ll_left = [l_fl] + left_edges + [-l_bl, -l_bbl]
    s_left = gmsh.model.geo.addPlaneSurface([gmsh.model.geo.addCurveLoop(ll_left)])

    # Create right face
    ll_right = [l_bbr, -l_fr] + right_edges + [l_br]
    s_right = gmsh.model.geo.addPlaneSurface([gmsh.model.geo.addCurveLoop(ll_right)])

    # Create substrate volume
    list_of_meta_surf = ( 
        black + white + corner_disks 
        + [item for sublist in x_edge_disks for item in sublist] 
        +  [item for sublist in y_edge_disks for item in sublist] 
        + [item for sublist in inner_disks for item in sublist] 
    ) 
    sl = gmsh.model.geo.addSurfaceLoop([s_bot, s_front, s_back, s_left, s_right] + list_of_meta_surf)
    subst_vol= gmsh.model.geo.addVolume([sl])

    # Create the source plane and first vacuum region
    # ***********************************************
    z_min = 0
    z_max = source_distance

    # Create the source plane
    p_vac1_fl = gmsh.model.geo.addPoint(x_min, y_min, z_max, hs)
    p_vac1_fr = gmsh.model.geo.addPoint(x_max, y_min, z_max, hs)
    p_vac1_br = gmsh.model.geo.addPoint(x_max, y_max, z_max, hs)
    p_vac1_bl = gmsh.model.geo.addPoint(x_min, y_max, z_max, hs)

    l_vac1_f = gmsh.model.geo.addLine(p_vac1_fl, p_vac1_fr)
    l_vac1_r = gmsh.model.geo.addLine(p_vac1_fr, p_vac1_br)
    l_vac1_b = gmsh.model.geo.addLine(p_vac1_br, p_vac1_bl)
    l_vac1_l = gmsh.model.geo.addLine(p_vac1_bl, p_vac1_fl)

    ll_vac1 = gmsh.model.geo.addCurveLoop([l_vac1_f, l_vac1_r, l_vac1_b, l_vac1_l])
    s_source = gmsh.model.geo.addPlaneSurface([ll_vac1])

    # Vertical edges connecting source plane to metasurface corners
    v1_fl = gmsh.model.geo.addLine(p_vac1_fl, corner_points[0])
    v1_fr = gmsh.model.geo.addLine(corner_points[1], p_vac1_fr)
    v1_bl = gmsh.model.geo.addLine(p_vac1_bl, corner_points[2])
    v1_br = gmsh.model.geo.addLine(corner_points[3], p_vac1_br)

    # Front face analogous to substrate front face
    ll_vac1_front = [-l_vac1_f, v1_fl] + front_edges + [v1_fr]
    s_vac1_front = gmsh.model.geo.addPlaneSurface([gmsh.model.geo.addCurveLoop(ll_vac1_front)])

    # Back face analogous to substrate back face
    ll_vac1_back = [l_vac1_b, v1_bl] + back_edges + [v1_br]
    s_vac1_back = gmsh.model.geo.addPlaneSurface([gmsh.model.geo.addCurveLoop(ll_vac1_back)])

    # Left face analogous to substrate left face
    ll_vac1_left = [v1_fl] + left_edges + [-v1_bl, l_vac1_l]
    s_vac1_left = gmsh.model.geo.addPlaneSurface([gmsh.model.geo.addCurveLoop(ll_vac1_left)])

    # Right face analogous to substrate right face
    ll_vac1_right = [-l_vac1_r, -v1_fr] + right_edges + [v1_br]
    s_vac1_right = gmsh.model.geo.addPlaneSurface([gmsh.model.geo.addCurveLoop(ll_vac1_right)])

    # Volume for first vacuum region: bounded below by metasurface (list_of_meta_surf) and above by source plane
    sl_vac1 = gmsh.model.geo.addSurfaceLoop([s_source, s_vac1_front, s_vac1_back, s_vac1_left, s_vac1_right] + list_of_meta_surf)
    vac1_vol = gmsh.model.geo.addVolume([sl_vac1])

    # Create the absorber plane and second vacuum region
    # **************************************************
    z_min = 0
    z_max = absorber_distance

    # Create the absorber plane
    p_vac2_fl = gmsh.model.geo.addPoint(x_min, y_min, z_max, ha)
    p_vac2_fr = gmsh.model.geo.addPoint(x_max, y_min, z_max, ha)
    p_vac2_br = gmsh.model.geo.addPoint(x_max, y_max, z_max, ha)
    p_vac2_bl = gmsh.model.geo.addPoint(x_min, y_max, z_max, ha)

    l_vac2_f = gmsh.model.geo.addLine(p_vac2_fl, p_vac2_fr)
    l_vac2_r = gmsh.model.geo.addLine(p_vac2_fr, p_vac2_br)
    l_vac2_b = gmsh.model.geo.addLine(p_vac2_br, p_vac2_bl)
    l_vac2_l = gmsh.model.geo.addLine(p_vac2_bl, p_vac2_fl)

    ll_vac2 = gmsh.model.geo.addCurveLoop([l_vac2_f, l_vac2_r, l_vac2_b, l_vac2_l])
    s_absorber = gmsh.model.geo.addPlaneSurface([ll_vac2])

    # Vertical edges connecting source plane to absorber plane
    v2_fl = gmsh.model.geo.addLine(p_vac1_fl, p_vac2_fl)
    v2_fr = gmsh.model.geo.addLine(p_vac1_fr, p_vac2_fr)
    v2_bl = gmsh.model.geo.addLine(p_vac1_bl, p_vac2_bl)
    v2_br = gmsh.model.geo.addLine(p_vac1_br, p_vac2_br)

    # Front face analogous to substrate front face
    ll_vac2_front = gmsh.model.geo.addCurveLoop([l_vac1_f, v2_fr, -l_vac2_f, -v2_fl])
    s_vac2_front = gmsh.model.geo.addPlaneSurface([ll_vac2_front])

    # Back face analogous to substrate back face
    ll_vac2_back = [-l_vac1_b, v2_br, l_vac2_b, -v2_bl]
    s_vac2_back = gmsh.model.geo.addPlaneSurface([gmsh.model.geo.addCurveLoop(ll_vac2_back)])

    # Left face analogous to substrate left face
    ll_vac2_left = [l_vac1_l, v2_fl, -l_vac2_l, -v2_bl]
    s_vac2_left = gmsh.model.geo.addPlaneSurface([gmsh.model.geo.addCurveLoop(ll_vac2_left)])   

    # Right face analogous to substrate right face
    ll_vac2_right = [-l_vac1_r, v2_fr, l_vac2_r, -v2_br]
    s_vac2_right = gmsh.model.geo.addPlaneSurface([gmsh.model.geo.addCurveLoop(ll_vac2_right)])

    # Volume for second vacuum region: bounded below by source plane and above by absorber plane
    sl_vac2 = gmsh.model.geo.addSurfaceLoop([s_absorber, s_vac2_front, s_vac2_back, s_vac2_left, s_vac2_right, s_source])
    vac2_vol = gmsh.model.geo.addVolume([sl_vac2])     

    # Synchronize the CAD kernel with the Gmsh model
    gmsh.model.geo.synchronize()


    # Periodic constraints
    # ********************
    # Define translation transformations for periodicity constraints
    translation_x = [
        1, 0, 0, 2*nx*wx,
        0, 1, 0, 0,
        0, 0, 1, 0,
        0, 0, 0, 1
    ]

    translation_y = [
        1, 0, 0, 0,
        0, 1, 0, 2*ny*wy,
        0, 0, 1, 0,
        0, 0, 0, 1
    ]

    # Set periodic constraints for substrate lateral faces
    gmsh.model.mesh.setPeriodic(2, [s_right], [s_left], translation_x)
    gmsh.model.mesh.setPeriodic(2, [s_back], [s_front], translation_y)

    # Set periodic constraints for vacuum region 1 lateral faces
    gmsh.model.mesh.setPeriodic(2, [s_vac1_right], [s_vac1_left], translation_x)
    gmsh.model.mesh.setPeriodic(2, [s_vac1_back], [s_vac1_front], translation_y)

    # Set periodic constraints for vacuum region 2 lateral faces
    gmsh.model.mesh.setPeriodic(2, [s_vac2_right], [s_vac2_left], translation_x)
    gmsh.model.mesh.setPeriodic(2, [s_vac2_back], [s_vac2_front], translation_y)


    # Return lists of geometric entities for physical group assignments
    # *****************************************************************
    substrate = [s_bot, s_front, s_back, s_left, s_right, subst_vol]
    metasurface = [black, white, corner_disks, x_edge_disks, y_edge_disks, inner_disks]

    vacuum1 = [s_source, s_vac1_front, s_vac1_back, s_vac1_left, s_vac1_right, vac1_vol]
    vacuum2 = [s_absorber, s_vac2_front, s_vac2_back, s_vac2_left, s_vac2_right, vac2_vol]

    return substrate, metasurface, vacuum1, vacuum2


def main():
    # Initialize Gmsh and create the model
    gmsh.initialize()
    gmsh.model.add(model_name)

    # Generate the model via geometric bottom-up construction
    substrate, metasurface, vacuum1, vacuum2 = meta_model(
        rect_x, rect_y, circle_r, 
        num_x, num_y,
        substrate_thickness,
        source_distance,
        absorber_distance,
        hx, hy, hc, hb, hs, ha 
    )

    # Define physical groups for the metasurface components
    gmsh.model.addPhysicalGroup(2, metasurface[0], tag=1, name="black_rectangles")
    gmsh.model.addPhysicalGroup(2, metasurface[1], tag=2, name="white_rectangles")
    gmsh.model.addPhysicalGroup(2, metasurface[2], tag=10, name="corner_disks")
    edge_disks = []
    for i, cut in enumerate(metasurface[3]): 
        # gmsh.model.addPhysicalGroup(2, cut, tag=100+i, name=f"x_edge_disk_{i}")
        edge_disks += cut
    for i, cut in enumerate(metasurface[4]):
        # gmsh.model.addPhysicalGroup(2, cut, tag=200+i, name=f"y_edge_disk_{i}")
        edge_disks += cut
    gmsh.model.addPhysicalGroup(2, edge_disks, tag=100, name="edge_disks")
    inner_disks = []
    for i, cut in enumerate(metasurface[5]):
        # gmsh.model.addPhysicalGroup(2, cut, tag=1000+i, name=f"inner_disk_{i}")
        inner_disks += cut
    gmsh.model.addPhysicalGroup(2, inner_disks, tag=1000, name="inner_disks")

    # Define physical groups for the substrate
    gmsh.model.addPhysicalGroup(2, [substrate[0]], tag=5000, name="substrate_bottom")
    gmsh.model.addPhysicalGroup(2, [substrate[1]], tag=5001, name="substrate_front")
    gmsh.model.addPhysicalGroup(2, [substrate[2]], tag=5002, name="substrate_back")
    gmsh.model.addPhysicalGroup(2, [substrate[3]], tag=5003, name="substrate_left")
    gmsh.model.addPhysicalGroup(2, [substrate[4]], tag=5004, name="substrate_right")
    gmsh.model.addPhysicalGroup(3, [substrate[5]], tag=6000, name="substrate_volume")

    # Define physical groups for the vacuum1 region
    gmsh.model.addPhysicalGroup(2, [vacuum1[0]], tag=20000, name="source_plane")
    gmsh.model.addPhysicalGroup(2, [vacuum1[1]], tag=20001, name="vacuum1_front")
    gmsh.model.addPhysicalGroup(2, [vacuum1[2]], tag=20002, name="vacuum1_back")
    gmsh.model.addPhysicalGroup(2, [vacuum1[3]], tag=20003, name="vacuum1_left")
    gmsh.model.addPhysicalGroup(2, [vacuum1[4]], tag=20004, name="vacuum1_right")
    gmsh.model.addPhysicalGroup(3, [vacuum1[5]], tag=21000, name="vacuum1_volume")

    # Define physical groups for the vacuum2 region
    gmsh.model.addPhysicalGroup(2, [vacuum2[0]], tag=30000, name="absorber_plane")
    gmsh.model.addPhysicalGroup(2, [vacuum2[1]], tag=30001, name="vacuum2_front")
    gmsh.model.addPhysicalGroup(2, [vacuum2[2]], tag=30002, name="vacuum2_back")
    gmsh.model.addPhysicalGroup(2, [vacuum2[3]], tag=30003, name="vacuum2_left")
    gmsh.model.addPhysicalGroup(2, [vacuum2[4]], tag=30004, name="vacuum2_right")
    gmsh.model.addPhysicalGroup(3, [vacuum2[5]], tag=31000, name="vacuum2_volume")

    # Set 2D meshing algorithm - 6 is frontal-Delaunay
    gmsh.option.setNumber("Mesh.Algorithm", 6)
    # Set 3D meshing algorithm - 4 is frontal
    gmsh.option.setNumber("Mesh.Algorithm3D", 4)

    # Enable mesh optimization
    gmsh.option.setNumber("Mesh.Optimize", 1)
    # Specifically enable Netgen optimization x3  
    gmsh.option.setNumber("Mesh.OptimizeNetgen", 3)  

    # Generate the mesh
    gmsh.model.mesh.generate()

    # Set mesh format to MSH 2.2 (ASCII)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.Binary", 0)  # 0 = ASCII, 1 = Binary

    # Save the mesh to file (robust path handling: ~ expansion, relative paths, auto-create dir)
    if save_mesh:
        out_dir = os.path.abspath(os.path.expanduser(path_to_mesh))
        os.makedirs(out_dir, exist_ok=True)
        gmsh.write(os.path.join(out_dir, f"{model_name}.msh"))

    # Display the generated model if desired
    if display_model:
        gmsh.fltk.run()

    # Terminate Gmsh
    gmsh.finalize()

    # Return model data
    return [substrate, metasurface, vacuum1, vacuum2]

if __name__ == "__main__":
    model_data = main()