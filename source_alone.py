# Author: Thierry Valet
# Copyright (C) 2025 - All rights reserved

"""
This script creates a source plane in between two absorber planes.
We use the native Gmsh kernel for the geometry definition.
A 3D mesh is generated for full-wave electromagnetic simulations.
--- ALL DIMENSIONS IN MILLIMETERS --- SOURCE AT z=0, CENTERED AT (x,y) ORIGIN ---
--- ABSORBER ABOVE AND BELOW---
"""

# Import necessary modules
import sys
import os
import gmsh

# General parameters
# Save mesh directly into the Palace test case folder for convenience
path_to_mesh = '/home/ubuntu/palace/sandbox/metasurfaces/test1'  # absolute path to save the mesh file
display_model = False  # Set to False for headless runs (no GUI)
save_mesh = True  # Set to True to save the mesh to file

# Model parameters
model_name = 'source_alone_v2'  # name of the model and mesh file
size_x = 60.     # rectangle dimension along x
size_y = 60.     # rectangle dimension along y
h = 5.0          # mesh size at absorber plane corners
absorber_distance = 75.0    # distance from source to absorber plane


# Funtion to create the full 3D model
def simple_model(
        sx, # Width of the rectangles
        sy, # Height of the rectangles
        d, # Distance from source to absorber planes
        h 
    ):

    # Create the source plane
    # ************************************
    x_min, x_max = -sx/2., sx/2.
    y_min, y_max = -sy/2., sy/2.
    p_fl = gmsh.model.geo.addPoint(x_min, y_min, 0., h)
    p_fr = gmsh.model.geo.addPoint(x_max, y_min, 0., h)
    p_br = gmsh.model.geo.addPoint(x_max, y_max, 0., h)
    p_bl = gmsh.model.geo.addPoint(x_min, y_max, 0., h)

    l_f = gmsh.model.geo.addLine(p_fr, p_fl)
    l_l = gmsh.model.geo.addLine(p_fl, p_bl)
    l_b = gmsh.model.geo.addLine(p_bl, p_br)
    l_r = gmsh.model.geo.addLine(p_br, p_fr)

    ll = gmsh.model.geo.addCurveLoop([l_f, l_l, l_b, l_r])
    s = gmsh.model.geo.addPlaneSurface([ll])

    # Create the bottom absorber plane
    # ************************************
    pb_fl = gmsh.model.geo.addPoint(x_min, y_min, -d, h)
    pb_fr = gmsh.model.geo.addPoint(x_max, y_min, -d, h)
    pb_br = gmsh.model.geo.addPoint(x_max, y_max, -d, h)
    pb_bl = gmsh.model.geo.addPoint(x_min, y_max, -d, h)

    lb_f = gmsh.model.geo.addLine(pb_fl, pb_fr)
    lb_r = gmsh.model.geo.addLine(pb_fr, pb_br)     
    lb_b = gmsh.model.geo.addLine(pb_br, pb_bl)
    lb_l = gmsh.model.geo.addLine(pb_bl, pb_fl) 

    ll_bot = gmsh.model.geo.addCurveLoop([lb_f, lb_r, lb_b, lb_l])
    s_bot = gmsh.model.geo.addPlaneSurface([ll_bot])

    # Create vacuum volume below the source plane
    # ************************************
    # Vertical edges connecting bottom absorber to source plane
    vb_fl = gmsh.model.geo.addLine(pb_fl, p_fl)
    vb_fr = gmsh.model.geo.addLine(pb_fr, p_fr)
    vb_bl = gmsh.model.geo.addLine(pb_bl, p_bl)    
    vb_br = gmsh.model.geo.addLine(pb_br, p_br)

    # Front face
    llb_vac_front = gmsh.model.geo.addCurveLoop([lb_f, vb_fr, l_f, -vb_fl])
    sb_vac_front = gmsh.model.geo.addPlaneSurface([llb_vac_front])

    # Back face
    llb_vac_back = gmsh.model.geo.addCurveLoop([-lb_b, vb_br, -l_b, -vb_bl])
    sb_vac_back = gmsh.model.geo.addPlaneSurface([llb_vac_back])

    # Left face
    llb_vac_left = gmsh.model.geo.addCurveLoop([-lb_l, vb_bl, -l_l, -vb_fl])
    sb_vac_left = gmsh.model.geo.addPlaneSurface([llb_vac_left])

    # Right face
    llb_vac_right = gmsh.model.geo.addCurveLoop([lb_r, vb_br, l_r, -vb_fr])
    sb_vac_right = gmsh.model.geo.addPlaneSurface([llb_vac_right])  

    # Volume
    vacb_sl = gmsh.model.geo.addSurfaceLoop([s_bot, sb_vac_front, sb_vac_back, sb_vac_left, sb_vac_right, s])
    vacb_vol = gmsh.model.geo.addVolume([vacb_sl])

    # Create the top absorber plane
    # ************************************
    pt_fl = gmsh.model.geo.addPoint(x_min, y_min, d, h)
    pt_fr = gmsh.model.geo.addPoint(x_max, y_min, d, h)
    pt_br = gmsh.model.geo.addPoint(x_max, y_max, d, h)     
    pt_bl = gmsh.model.geo.addPoint(x_min, y_max, d, h)

    lt_f = gmsh.model.geo.addLine(pt_fl, pt_fr)
    lt_r = gmsh.model.geo.addLine(pt_fr, pt_br)
    lt_b = gmsh.model.geo.addLine(pt_br, pt_bl)
    lt_l = gmsh.model.geo.addLine(pt_bl, pt_fl)

    ll_top = gmsh.model.geo.addCurveLoop([lt_f, lt_r, lt_b, lt_l])
    s_top = gmsh.model.geo.addPlaneSurface([ll_top])

    # Create vacuum volume above the source plane
    # ************************************
    # Vertical edges connecting top absorber to source plane
    vt_fl = gmsh.model.geo.addLine(p_fl, pt_fl)
    vt_fr = gmsh.model.geo.addLine(p_fr, pt_fr)
    vt_bl = gmsh.model.geo.addLine(p_bl, pt_bl)
    vt_br = gmsh.model.geo.addLine(p_br, pt_br)

    # Front face
    llt_vac_front = gmsh.model.geo.addCurveLoop([-lt_f, -vt_fl, -l_f, vt_fr])
    st_vac_front = gmsh.model.geo.addPlaneSurface([llt_vac_front])

    # Back face
    llt_vac_back = gmsh.model.geo.addCurveLoop([lt_b, -vt_bl, l_b, vt_br])
    st_vac_back = gmsh.model.geo.addPlaneSurface([llt_vac_back])

    # Left face
    llt_vac_left = gmsh.model.geo.addCurveLoop([lt_l, -vt_fl, l_l, vt_bl])
    st_vac_left = gmsh.model.geo.addPlaneSurface([llt_vac_left])

    # Right face
    llt_vac_right = gmsh.model.geo.addCurveLoop([-lt_r, -vt_fr, -l_r, vt_br])
    st_vac_right = gmsh.model.geo.addPlaneSurface([llt_vac_right])

    # Volume
    vact_sl = gmsh.model.geo.addSurfaceLoop([s_top, st_vac_front, st_vac_back, st_vac_left, st_vac_right, s])
    vact_vol = gmsh.model.geo.addVolume([vact_sl])


    # Synchronize the CAD kernel with the Gmsh model
    gmsh.model.geo.synchronize()


    # Periodic constraints
    # ********************
    # Define translation transformations for periodicity constraints
    translation_x = [
        1, 0, 0, sx,
        0, 1, 0, 0,
        0, 0, 1, 0,
        0, 0, 0, 1
    ]

    translation_y = [
        1, 0, 0, 0,
        0, 1, 0, sy,
        0, 0, 1, 0,
        0, 0, 0, 1
    ]

    # Set periodic constraints for bottom vacuum region lateral faces
    gmsh.model.mesh.setPeriodic(2, [sb_vac_right], [sb_vac_left], translation_x)
    gmsh.model.mesh.setPeriodic(2, [sb_vac_back], [sb_vac_front], translation_y)


    # Set periodic constraints for top vacuum region lateral faces
    gmsh.model.mesh.setPeriodic(2, [st_vac_right], [st_vac_left], translation_x)
    gmsh.model.mesh.setPeriodic(2, [st_vac_back], [st_vac_front], translation_y)


    # Return lists of geometric entities for physical group assignments
    # *****************************************************************
    front = [sb_vac_front, st_vac_front]
    back = [sb_vac_back, st_vac_back]
    left = [sb_vac_left, st_vac_left]
    right = [sb_vac_right, st_vac_right]

    return s, s_bot, s_top, front, back, left, right, vacb_vol, vact_vol


def main():
    # Initialize Gmsh and create the model
    gmsh.initialize()
    gmsh.model.add(model_name)

    # Generate the model via geometric bottom-up construction
    ( source, bottom_absorber, top_absorber,
    front, back, left, right, 
    bottom_vac, top_vac ) = simple_model(
        size_x, size_y,
        absorber_distance,
        h
    ) 

    # Synchronize the CAD kernel with the Gmsh model
    gmsh.model.geo.synchronize()

    # Define physical groups for periodic boundary faces
    gmsh.model.addPhysicalGroup(2, front, tag=101, name="periodic_front")
    gmsh.model.addPhysicalGroup(2, back, tag=102, name="periodic_back")
    gmsh.model.addPhysicalGroup(2, left, tag=103, name="periodic_left")
    gmsh.model.addPhysicalGroup(2, right, tag=104, name="periodic_right")

    # Define physical groups for the source
    gmsh.model.addPhysicalGroup(2, [source], tag=1000, name="source_plane")

    # Define physical groups for the absorbers
    gmsh.model.addPhysicalGroup(2, [bottom_absorber], tag=2000, name="bottom_absorber")
    gmsh.model.addPhysicalGroup(2, [top_absorber], tag=2001, name="top_absorber")

    # Define physical groups for the vacuum regions
    gmsh.model.addPhysicalGroup(3, [bottom_vac], tag=10000, name="bottom_vacuum")
    gmsh.model.addPhysicalGroup(3, [top_vac], tag=11000, name="top_vacuum")

    # Synchronize the CAD kernel with the Gmsh model
    gmsh.model.geo.synchronize()

    # Set 2D meshing algorithm - 6 is frontal-Delaunay
    gmsh.option.setNumber("Mesh.Algorithm", 6)
    # Set 3D meshing algorithm - 4 is frontal
    gmsh.option.setNumber("Mesh.Algorithm3D", 4)

    # Mesh optimization options 
    gmsh.option.setNumber("Mesh.Optimize", 0)
    gmsh.option.setNumber("Mesh.OptimizeNetgen", 0)

    # Generate the mesh
    gmsh.model.mesh.generate()

    # Set mesh format to MSH 2.2 (ASCII) to ensure compatibility with Palace/MFEM
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
    return 0

if __name__ == "__main__":
    model_data = main()