# Author: Thierry Valet
# Copyright (C) 2025 - All rights reserved

"""
This script creates a chequerboard metasurface unit cell,
and then a supercell by replicating the unit cell in x and y directions.
We use the native Gmsh kernel.
"""

# Import necessary modules
import sys
import numpy as np
import gmsh

# Function to create a sub-unit cell with one rectangle and 
# arc-circle cutouts at each corner
def sub_unit_cell(
        x0, # Bottom left corner x-coordinate
        y0, # Bottom left corner y-coordinate
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

    p4 = gmsh.model.geo.addPoint(x0+wx, y0, 0, hc)
    p5 = gmsh.model.geo.addPoint(x0+wx-r, y0, 0, hc)
    p6 = gmsh.model.geo.addPoint(x0+wx, y0+r, 0, hc)  

    l2 = gmsh.model.geo.addLine(p5, p4)
    c1 = gmsh.model.geo.addCircleArc(p6, p4, p5)
    l3 = gmsh.model.geo.addLine(p4, p6)
    ll1 = gmsh.model.geo.addCurveLoop([l2, l3, c1])
    s1 = gmsh.model.geo.addPlaneSurface([ll1])

    p7 = gmsh.model.geo.addPoint(x0+wx, y0+wy/2, 0, hy)

    p8 = gmsh.model.geo.addPoint(x0+wx, y0+wy, 0, hc)
    p9 = gmsh.model.geo.addPoint(x0+wx, y0+wy-r, 0, hc)
    p10 = gmsh.model.geo.addPoint(x0+wx-r, y0+wy, 0, hc)

    l4 = gmsh.model.geo.addLine(p9, p8)
    c2 = gmsh.model.geo.addCircleArc(p10, p8, p9)
    l5 = gmsh.model.geo.addLine(p8, p10)
    ll2 = gmsh.model.geo.addCurveLoop([l5, c2, l4])
    s2 = gmsh.model.geo.addPlaneSurface([ll2])

    p11 = gmsh.model.geo.addPoint(x0+wx/2, y0+wy, 0, hx)

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
    bottom = [l0, l8, l9, l2]
    right = [l3, l10, l11, l4]
    top = [l6, -l13, -l12, -l5]
    left = [-l1, -l15, -l14, -l7]

    return [surfaces , bottom, right, top, left]

# Function to create a sub-unit cell with one rectangle and 
# arc-circle cutouts at each corner
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
    c = [cb, cm, ct]

    # arrays of edges in the unit cell
    bottom = np.array(lower_left[1] + lower_right[1], dtype=int)
    right = np.array(lower_right[2] + upper_right[2], dtype=int)
    top = np.array(upper_left[3] + upper_right[3], dtype=int)
    left = np.array(lower_left[4] + upper_left[4], dtype=int)

    return [black_rectangles, white_rectangles, c, bottom, right, top, left]

def period(
        wx, # Width of the rectangles
        wy, # Height of the rectangles
        r,  # Radius of the circular regions
        hx, # Mesh size at rectangle x edges
        hy, # Mesh size at rectangle y edges
        hc,   # Mesh size at circular perimeter
        nx, # replication count in x
        ny  # replication count in y
    ):
    # Iterate over the number of replications in x and y
    black = []
    white = []
    for i in range(nx):
        for j in range(ny):
            # Create the unit cell at the appropriate position
            out = unit_cell((2*i - nx) * wx, (2*j - ny) * wy, wx, wy, r, hx, hy, hc)
            black += out[0]
            white += out[1]


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
substrate_thickness = 0.01
source_distance = 0.06
absorber_distance = 0.06

# Initialize Gmsh and create the model
gmsh.initialize()
gmsh.model.add(model_name)

# Generate the periodic structure
period(rect_x, rect_y, circle_r, hx, hy, hc, num_x, num_y)

# Remove duplicate entities
gmsh.model.geo.removeAllDuplicates()

gmsh.model.geo.synchronize()

# Display the generated metasurface
gmsh.fltk.run()

gmsh.finalize()