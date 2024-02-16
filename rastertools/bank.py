import numpy as np
import rioxarray
import geopandas as gpd
from shapely.geometry import LineString
from .boundary import *
from .centerline import *


def delineate_bank_lines(dem_file: str) -> gpd.GeoDataFrame:
    dem = rioxarray.open_rasterio(dem_file, masked=True)
    utm_crs = dem.rio.estimate_utm_crs()
    dem = dem.rio.reproject(utm_crs)

    ## determine banks by the slopes and percentile of elevation
    dx = dem.differentiate('x')
    dy = dem.differentiate('y')
    x_steep = np.nanpercentile(np.abs(dx.values.flatten()), 95)
    y_steep = np.nanpercentile(np.abs(dy.values.flatten()), 95)

    boundary = dem.where((np.abs(dx) > x_steep) & (np.abs(dy) > y_steep), drop=True)
    dem_lb = np.nanpercentile(boundary, 2.5)
    dem_ub = np.nanpercentile(boundary, 5)
    banks = dem.where((dem > dem_lb) & (dem < dem_ub), drop=True)
    r = banks.rolling(x=3, y=3, center=True, min_periods=1)
    r_max = r.max(skipna=True)
    thin_banks = banks.where(banks == r_max, drop=True)
    
    ## convert banks from raster to a vector
    thin_banks.name = 'elevation'
    thin_banks = thin_banks.to_dataframe()
    thin_banks.dropna(inplace=True)
    thin_banks.reset_index(inplace=True)
    thin_banks = thin_banks[['x', 'y', 'elevation']]
    bd_points = gpd.GeoDataFrame(thin_banks, geometry=gpd.points_from_xy(thin_banks.x, thin_banks.y), crs=utm_crs)

    ## dem -> boundary polygon -> centerline
    boundary = delineate_boundary(dem_file)
    boundary_geom = boundary.geometry[0]
    centerline = delineate_centerline(boundary_geom)

    ## get SN coordinates to sort and split the bank points
    cord_converter = CoordConverter(centerline)
    getSN = lambda row: cord_converter.xy2sn_coord(row['x'], row['y'])
    bd_points['s'], bd_points['n'], _ = zip(*bd_points.apply(getSN, axis=1))
    bd_points.sort_values('s', inplace=True)
    left_bank_points = bd_points[bd_points['n'] < 0]
    right_bank_points = bd_points[bd_points['n'] > 0]

    ## form two bank lines and store in a GeoDataFrame
    left_bank_line = LineString(left_bank_points[['x', 'y']].values)
    right_bank_line = LineString(right_bank_points[['x', 'y']].values)
    bank_lines = gpd.GeoDataFrame(geometry=gpd.GeoSeries([left_bank_line, right_bank_line]), crs=utm_crs)

    return bank_lines

