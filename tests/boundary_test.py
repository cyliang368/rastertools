#%%
from rastertools import *
import os
from glob import glob
import warnings
warnings.filterwarnings("ignore")

dem_files = glob("../data/*/input/dem_clip.tif")

for dem_file in dem_files:
    bd_gdf = delineate_boundary(dem_file=dem_file, threshold=2.5)
    raster = rioxarray.open_rasterio(dem_file, masked=True)
    utm_crs = raster.rio.estimate_utm_crs()
    raster = raster.rio.reproject(utm_crs)
    
    fig, ax = plt.subplots()
    raster.plot(ax=ax)
    bd_gdf.plot(ax=ax, facecolor='none', edgecolor='r')
    plt.tight_layout()
    plt.show()
# %%
