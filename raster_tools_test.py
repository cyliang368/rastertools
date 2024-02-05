#%%
import xarray as xr
import rioxarray
import geopandas as gpd
import pandas as pd
import matplotlib.pyplot as plt
import os
from shapely.geometry import *
from osgeo import gdal, osr

bd_file = '../XS_processing/xs_tools/Brazos/input/boundary.shp'
mesh_points = '../XS_processing/xs_tools/Brazos/output/shp/mesh_points.shp'
bd_gdf = gpd.read_file(bd_file)
mesh_gdf = gpd.read_file(mesh_points)

mesh_surface = './temp/mesh_surface.tif'
mesh_gdf = gpd.read_file(mesh_points)
xmin, ymin, xmax, ymax = mesh_gdf.loc[0, 'geometry'].bounds
resolution = 5
w = abs(xmax - xmin) // resolution
h = abs(ymax - ymin) // resolution
rasterDs = gdal.Grid(mesh_surface, mesh_points, format='GTiff', algorithm='invdist', zfield='z', width=w, height=h)
rasterDs = None

#%%
mesh_gdf.total_bounds

#%%
from rasterio.enums import Resampling

DEM_file = "./temp/USGS_13_n30w096_20211103_clip_proj.tif"

DEM = rioxarray.open_rasterio(DEM_file, masked=True)
xRes, yRes = DEM.rio.resolution()
scale_factor = xRes / resolution 
new_width = int(DEM.rio.width * scale_factor)
new_height = int(DEM.rio.height * scale_factor)
DEM = DEM.rio.reproject(
    DEM.rio.crs,
    shape=(new_height, new_width),
    resampling=Resampling.bilinear,
)
DEM.plot()

#%%
xds2 = rioxarray.open_rasterio("./temp/mesh_surface.tif", masked=True)
rest_DEM = xds2.rio.clip(bd_gdf.geometry, bd_gdf.crs, drop=True)
mesh_clip_file = './temp/mesh_clip.tif'
rest_DEM.rio.to_raster(mesh_clip_file)
rest_DEM.plot()

#%%
import numpy as np
reproject_mesh =rest_DEM.rio.reproject_match(DEM, nodata=np.nan)
reproject_mesh

merged_DEM = DEM.where(reproject_mesh.isnull(), reproject_mesh)
merged_DEM.rio.to_raster('./temp/merged.tif')
merged_DEM.plot() 



# %%
