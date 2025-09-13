# Author: Thierry Valet
# Copyright (C) 2025 - All rights reserved

"""
This script creates a chequerboard metasurface unit cell and a supercell by replicating the unit cell in x and y directions.
"""
# Import necessary modules
import gmsh, sys

# Function to create a unit cell with rectangles and circular cutouts
def unit_cell(
        x0, # Bottom left corner x-coordinate
        y0, # Bottom left corner y-coordinate
        wx, # Width of the rectangles
        wy, # Height of the rectangles
        r,  # Radius of the circular region
        h   # Mesh size at points
    ):
    # Creates lower left rectangle
    rect0 = gmsh.model.occ.addRectangle(x0, y0, 0, wx, wy)
    # Creates lower right rectangle
    rect1 = gmsh.model.occ.addRectangle(x0+wx, y0, 0, wx, wy)
    # Creates upper left rectangle
    rect2 = gmsh.model.occ.addRectangle(x0, y0+wy, 0, wx, wy)
    # Creates upper right rectangle
    rect3 = gmsh.model.occ.addRectangle(x0+wx, y0+wy, 0, wx, wy)


    # Creates a disk at lower left corner of the rectangle
    disk0 = gmsh.model.occ.addDisk(x0, y0, 0, r, r)
    # Compute the intersection of the rectangle and the disk and remove the disk
    dimTag0, dimTagMap = gmsh.model.occ.intersect([(2,rect0)], [(2,disk0)], removeObject=False)
    # Cut the rectangle with the intersection
    rect0_cut, dimTagMap = gmsh.model.occ.cut([(2,rect0)], dimTag0, removeObject=True, removeTool=False)

    # Creates a disk at the lower connecting corner of the rectangles
    disk1 = gmsh.model.occ.addDisk(x0+wx, y0, 0, r, r)
    # Compute the intersection of the rectangles and the disk and remove the disk
    dimTag1, dimTagMap = gmsh.model.occ.intersect(rect0_cut + [(2,rect1)], [(2,disk1)], removeObject=False)
    # Cut the rectangles with the intersection
    rect0_cut, dimTagMap = gmsh.model.occ.cut(rect0_cut, dimTag1, removeObject=True, removeTool=False)
    rect1_cut, dimTagMap = gmsh.model.occ.cut([(2,rect1)], dimTag1, removeObject=True, removeTool=False)

    # Creates a disk at the lower right corner of the rectangles
    disk2 = gmsh.model.occ.addDisk(x0+2*wx, y0, 0, r, r)
    # Compute the intersection of the rectangle and the disk and remove the disk
    dimTag2, dimTagMap = gmsh.model.occ.intersect(rect1_cut, [(2,disk2)], removeObject=False)
    # Cut the rectangle with the intersection
    rect1_cut, dimTagMap = gmsh.model.occ.cut(rect1_cut, dimTag2, removeObject=True, removeTool=False)    

    # Creates a disk at the left connecting corner of the rectangles
    disk3 = gmsh.model.occ.addDisk(x0, y0+wy, 0, r, r)
    # Compute the intersection of the rectangles and the disk and remove the disk
    dimTag3, dimTagMap = gmsh.model.occ.intersect(rect0_cut+[(2,rect2)], [(2,disk3)], removeObject=False)
    # Cut the rectangles with the intersection
    rect0_cut, dimTagMap = gmsh.model.occ.cut(rect0_cut, dimTag3, removeObject=True, removeTool=False)
    rect2_cut, dimTagMap = gmsh.model.occ.cut([(2,rect2)], dimTag3, removeObject=True, removeTool=False)

    # Creates a disk at the center connecting corner of the rectangles
    disk4 = gmsh.model.occ.addDisk(x0+wx, y0+wy, 0, r, r)
    # Compute the intersection of the rectangles and the disk and remove the disk
    dimTag4, dimTagMap = gmsh.model.occ.intersect(rect0_cut + rect1_cut + rect2_cut + [(2,rect3)], [(2,disk4)], removeObject=False)
    # Cut the rectangles with the intersection
    rect0_cut, dimTagMap = gmsh.model.occ.cut(rect0_cut, dimTag4, removeObject=True, removeTool=False)
    rect1_cut, dimTagMap = gmsh.model.occ.cut(rect1_cut, dimTag4, removeObject=True, removeTool=False)
    rect2_cut, dimTagMap = gmsh.model.occ.cut(rect2_cut, dimTag4, removeObject=True, removeTool=False)
    rect3_cut, dimTagMap = gmsh.model.occ.cut([(2,rect3)], dimTag4, removeObject=True, removeTool=False)

    # Creates a disk at the right connecting corner of the rectangles
    disk5 = gmsh.model.occ.addDisk(x0+2*wx, y0+wy, 0, r, r)
    # Compute the intersection of the rectangles and the disk and remove the disk
    dimTag5, dimTagMap = gmsh.model.occ.intersect(rect1_cut + rect3_cut, [(2,disk5)], removeObject=False)
    # Cut the rectangles with the intersection
    rect1_cut, dimTagMap = gmsh.model.occ.cut(rect1_cut, dimTag5, removeObject=True, removeTool=False)
    rect3_cut, dimTagMap = gmsh.model.occ.cut([(2,rect3)], dimTag5, removeObject=True, removeTool=False)

    # Creates a disk at the upper left corner of the rectangles
    disk6 = gmsh.model.occ.addDisk(x0, y0+2*wy, 0, r, r)
    # Compute the intersection of the rectangle and the disk and remove the disk
    dimTag6, dimTagMap = gmsh.model.occ.intersect(rect2_cut, [(2,disk6)], removeObject=False)
    # Cut the rectangle with the intersection
    rect2_cut, dimTagMap = gmsh.model.occ.cut([(2,rect2)], dimTag6, removeObject=True, removeTool=False)

    # Creates a disk at the upper connecting corner of the rectangles
    disk7 = gmsh.model.occ.addDisk(x0+wx, y0+2*wy, 0, r, r)
    # Compute the intersection of the rectangles and the disk and remove the disk
    dimTag7, dimTagMap = gmsh.model.occ.intersect(rect2_cut + rect3_cut, [(2,disk7)], removeObject=False)
    # Cut the rectangles with the intersection
    rect2_cut, dimTagMap = gmsh.model.occ.cut(rect2_cut, dimTag7, removeObject=True, removeTool=False)
    rect3_cut, dimTagMap = gmsh.model.occ.cut(rect3_cut, dimTag7, removeObject=True, removeTool=False)

    # Creates a disk at the upper right corner of the rectangles
    disk8 = gmsh.model.occ.addDisk(x0+2*wx, y0+2*wy, 0, r, r)
    # Compute the intersection of the rectangle and the disk and remove the disk
    dimTag8, dimTagMap = gmsh.model.occ.intersect(rect3_cut, [(2,disk8)], removeObject=False)
    # Cut the rectangle with the intersection
    rect3_cut, dimTagMap = gmsh.model.occ.cut(rect3_cut, dimTag8, removeObject=True, removeTool=False)

    return [rect0_cut, rect1_cut, rect2_cut, rect3_cut], [dimTag0, dimTag1, dimTag2, dimTag3, dimTag4, dimTag5, dimTag6, dimTag7, dimTag8]

# Function to create a supercell by replicating the unit cell
def supercell(
        nx, # Number of unit cells in x-direction
        ny, # Number of unit cells in y-direction
        wx, # Width of the unit cell
        wy, # Height of the unit cell
        r,  # Radius of the disks
        h   # Mesh size at points
        ):
    list_of_rectangles0 = []
    list_of_rectangles1 = []
    list_of_disks = []
    for i in range(nx):
        for j in range(ny):
            x0 = (2*i - nx)*wx
            y0 = (2*j-ny)*wy
            rectangles, disks = unit_cell(x0, y0, wx, wy, r, h)
            list_of_rectangles0 += [rectangles[0][0][1]]
            list_of_rectangles1 += [rectangles[1][0][1]]
            list_of_rectangles1 += [rectangles[2][0][1]]
            list_of_rectangles0 += [rectangles[3][0][1]]
            list_of_disks += [disks[0][0][1]]
            list_of_disks += [disks[1][0][1]]
            list_of_disks += [disks[1][1][1]]
            list_of_disks += [disks[2][0][1]]
            list_of_disks += [disks[3][0][1]]
            list_of_disks += [disks[3][1][1]]
            list_of_disks += [disks[4][0][1]]
            list_of_disks += [disks[4][1][1]]
            list_of_disks += [disks[4][2][1]]
            list_of_disks += [disks[4][3][1]]
            list_of_disks += [disks[5][0][1]]
            list_of_disks += [disks[5][1][1]]
            list_of_disks += [disks[6][0][1]]
            list_of_disks += [disks[7][0][1]]
            list_of_disks += [disks[7][1][1]]
            list_of_disks += [disks[8][0][1]]
    return list_of_rectangles0, list_of_rectangles1, list_of_disks

# Configuration parameters
model_name = 'meta_surf'
num_x = 5          # replication count in x
num_y = 5          # replication count in y
rect_x = 0.003      # rectangle dimension along x
rect_y = 0.003      # rectangle dimension along y
circle_r = 0.0005   # disk radius
mesh_h0 = 1e-4
substrate_thickness = 0.01
source_distance = 0.06
absorber_distance = 0.06

# Initialize Gmsh and create the model
gmsh.initialize()
gmsh.model.add(model_name)

# Generate the supercell
list_of_rectangles_b, list_of_rectangles_w, list_of_disks = supercell(num_x, num_y, rect_x, rect_y, circle_r, mesh_h0)

# Generate the substrate as a box
substrate = gmsh.model.occ.addBox(-num_x*rect_x, -num_y*rect_y, -substrate_thickness, 2*num_x*rect_x, 2*num_y*rect_y, substrate_thickness)

# Generate the first vacuum layer above the metasurface and surface current plane
vacuum_layer0 = gmsh.model.occ.addBox(-num_x*rect_x, -num_y*rect_y, 0, 2*num_x*rect_x, 2*num_y*rect_y, source_distance)

# Generate the second vacuum layer above the first vacuum layer and absorber plane
vacuum_layer1 = gmsh.model.occ.addBox(-num_x*rect_x, -num_y*rect_y, source_distance, 2*num_x*rect_x, 2*num_y*rect_y, absorber_distance)

# Remove duplicates
gmsh.model.occ.removeAllDuplicates()

# Synchronize the CAD kernel with the Gmsh model
gmsh.model.occ.synchronize()

# Define physical groups for the metasurface components
gmsh.model.addPhysicalGroup(2, list_of_rectangles_b, tag=1, name="rectangles_b")
gmsh.model.addPhysicalGroup(2, list_of_rectangles_w, tag=2, name="rectangles_w")
gmsh.model.addPhysicalGroup(2, list_of_disks, tag=3, name="disks")

# Set 2D meshing algorithm - 6 is frontal-Delaunay
gmsh.option.setNumber("Mesh.Algorithm", 6)
# Set 3D meshing algorithm - 4 is frontal
gmsh.option.setNumber("Mesh.Algorithm3D", 4)

# Generate the mesh
gmsh.model.mesh.generate()

# Display the generated metasurface
gmsh.fltk.run()

gmsh.finalize()