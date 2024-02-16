from __future__ import annotations

import os
import geopandas as gpd
import pandas as pd
import numpy as np
import rasterio as rio
import rioxarray
import xarray as xr
from shapely.ops import unary_union
import geopandas as gpd
from shapely.geometry import Polygon
import matplotlib.pyplot as plt
import warnings
from joblib import Parallel, delayed


def delineate_boundary(dem_file: str, threshold: float = 2.5) -> gpd.GeoDataFrame:
    """delineate boundary by DEM based on slopes and percentile threshold of elevation

    Parameters
    ----------
    dem_files : str
        _description_
    threshold : float, optional
        a percentile of elevation used to determine the water body a higher value will derive larger area, by default 2.5
    """

    ## read dem and reproject to UTM
    raster = rioxarray.open_rasterio(dem_file, masked=True)
    utm_crs = raster.rio.estimate_utm_crs()
    raster = raster.rio.reproject(utm_crs)

    ## determine the threshold of the water body by slopes and percentile of elevation
    dx = raster.differentiate('x')
    dy = raster.differentiate('y')
    x_steep = np.nanpercentile(np.abs(dx.values.flatten()), 95)
    y_steep = np.nanpercentile(np.abs(dy.values.flatten()), 95)
    boundary = raster.where((np.abs(dx) > x_steep) & (np.abs(dy) > y_steep), drop=True)
    dem_threshold = np.nanpercentile(boundary, threshold)

    ## convert the dem to a vector
    x, y, elevation = raster.x.values, raster.y.values, raster.values
    x, y = np.meshgrid(x, y)
    x, y, elevation = x.flatten(), y.flatten(), elevation.flatten()

    dem_pd = pd.DataFrame.from_dict({'elevation': elevation, 'x': x, 'y': y})
    dem_pd = dem_pd[dem_pd['elevation'] < dem_threshold]
    dem_vector = gpd.GeoDataFrame(geometry=gpd.GeoSeries.from_xy(dem_pd['x'], dem_pd['y'], crs=raster.rio.crs))
    dem_vector = dem_vector.buffer(5, cap_style=3)
    dem_vector = dem_vector.to_crs(utm_crs)
    geom_arr = []
    
    # converting GeoSeries to list of geometries
    geoms = list(dem_vector)

    # converting geometries list to nested list of geometries
    geom_arr = [geoms[i:i+10000] for i in range(0, len(geoms), 10000)]

    # perform union operation of chunks of geometries
    geom_union = Parallel(n_jobs=-1, verbose=0)(delayed(unary_union)(geom) for geom in geom_arr)
    
    # Perform union operation on returned unioned geometries
    total_union = unary_union(geom_union)
    union_vector_gdf = gpd.GeoDataFrame(geometry=gpd.GeoSeries(total_union))

    # fig, ax = plt.subplots()
    # raster.plot(ax=ax)
    # union_vector_gdf.buffer(15).plot(ax=ax, facecolor='none', edgecolor='r')
    # plt.show()
    return union_vector_gdf