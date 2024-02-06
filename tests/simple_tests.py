from rastertools import *
import os
from glob import glob
import warnings
warnings.filterwarnings("ignore")

test_folders = glob("../data/*")
for folder in test_folders:
    input_folder = folder + "/input"
    output_folder = folder + "/output"

    mesh_file = input_folder + "/mesh_points.geojson"
    boundary_file = input_folder + "/boundary.geojson"
    dem_file = input_folder + "/dem_clip.tif"
    surface_file = output_folder + "/surface.tif"
    merged_file = output_folder + "/merged.tif"
    
    
    value_type = 'depth-' if ('brazos' in folder) else 'depth+'

    interpolate_points_to_raster(mesh_file=mesh_file, 
                                 surface_file=surface_file, 
                                 boundary_file=boundary_file, 
                                 resolution=5)
    
    
    integrate_surface_with_dem(surface_file=surface_file,
                               dem_file=dem_file,
                               merged_file=merged_file, 
                               value_type=value_type, 
                               resolution=5)
        
    