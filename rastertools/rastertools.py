#%%
import geopandas as gpd
import numpy as np
from scipy import interpolate
import rasterio as rio
import rioxarray
from rasterio.transform import Affine
from rasterio.enums import Resampling

def clip_raster_by_mesh_extent(mesh_file, raster_file, output_path, dest_crs=None):
    ## check mesh file
    if (mesh_file is None):
        raise "Please set up a mesh_file before using clip_dem_by_input_extent!"
    elif (isinstance(mesh_file, str)):
        mesh_gdf = gpd.read_file(mesh_file)
    elif (isinstance(mesh_file, gpd.GeoDataFrame)):
        mesh_gdf = mesh_file
    
    ## check DEM file
    if (raster_file is None):
        raise "Please set up a dem_file before using clip_dem_by_input_extent!"
    elif (isinstance(raster_file, str)):
        raster = rioxarray.open_rasterio(raster_file, masked=True)
    elif (isinstance(raster_file, rioxarray.raster.RasterDataArray)):
        raster = raster_file

    ## get the extent of the mesh
    minx, miny, maxx, maxy = mesh_gdf.to_crs(raster.rio.crs).total_bounds

    ## clip the DEM by the extent of the mesh
    clipped = raster.rio.clip_box(minx, miny, maxx, maxy)

    ## reproject the clipped DEM to the desired CRS
    if (dest_crs is not None):
        clipped = clipped.rio.reproject(dest_crs)

    ## save the clipped DEM to the output path
    if (output_path is not None):
        clipped.rio.to_raster(output_path)
    
    return clipped

def reproject_raster(raster_file, dest_crs, output_file=None):
    raster = rioxarray.open_rasterio(raster_file, masked=True)
    raster = raster.rio.reproject(dest_crs)
    if (output_file is not None):
        raster.rio.to_raster(output_file)
    return raster

def interpolate_points_to_raster(mesh_file, surface_file=None, boundary_file=None, resolution=5):
    ## check mesh file
    if (mesh_file is None):
        raise "Please set up a mesh_file before using clip_dem_by_input_extent!"
    elif (isinstance(mesh_file, str)):
        mesh_gdf = gpd.read_file(mesh_file)
    elif (isinstance(mesh_file, gpd.GeoDataFrame)):
        mesh_gdf = mesh_file

    ## calculate the resolution of the interpolated raster
    minx, miny, maxx, maxy = mesh_gdf.total_bounds
    w = abs(maxx - minx) // resolution
    h = abs(maxy - miny) // resolution

    vx = np.arange(minx, maxx+resolution, resolution)
    vy = np.arange(miny, maxy+resolution, resolution)
    gridX, gridY = np.meshgrid(vx, vy)
    XYZ = mesh_gdf.geometry.apply(lambda geom: geom.coords[0])
    x, y, z = zip(*XYZ)
    points = list(zip(x,y))
    interp_grid = interpolate.griddata(points, z, (gridX, gridY), method='cubic')

    transform = Affine.translation(gridX[0][0]-resolution/2, gridY[0][0]-resolution/2)*Affine.scale(resolution,resolution)

    with rio.open(
            surface_file, 'w',
            driver='GTiff',
            dtype=rio.float32,
            count=1,
            width=w,
            height=h,
            crs=mesh_gdf.crs,
            transform = transform) as dst:
        dst.write(interp_grid, indexes=1)

    ## clip the raster by the boundary or the extent of the mesh
    surface = rioxarray.open_rasterio(surface_file, masked=True)
    if (boundary_file is not None):
        bd_gdf = gpd.read_file(boundary_file)
        bd_gdf.to_crs(surface.rio.crs)
        surface = surface.rio.clip(bd_gdf.geometry, bd_gdf.crs, drop=True)
    else:
        surface = surface.rio.clip_box(minx, miny, maxx, maxy)

    if (surface_file is not None):
        surface.rio.to_raster(surface_file)
    return surface          

def integrate_surface_with_dem(surface_file, dem_file, merged_file, value_type='elevation', resolution=None):
    ## open two rasters that need to be integrated
    dem = rioxarray.open_rasterio(dem_file, masked=True)
    surface = rioxarray.open_rasterio(surface_file, masked=True)

    ## calculate the target resolution
    if (resolution is not None):
        ## reproject to utm crs to make sure the unit is in meters
        utm_crs = surface.rio.estimate_utm_crs()
        surface = surface.rio.reproject(utm_crs)
        dem = dem.rio.reproject(utm_crs)

        ## calculate the scale factor by the original and target resolutions
        xRes, yRes = dem.rio.resolution()
        scale_factor = xRes / resolution
        new_width = int(dem.rio.width * scale_factor)
        new_height = int(dem.rio.height * scale_factor)
        dem = dem.rio.reproject(
            surface.rio.crs,
            shape = (new_height, new_width),
            resampling = Resampling.bilinear,
        )
    else:
        dem = dem.rio.reproject_match(surface)

    reproject_surface = surface.rio.reproject_match(dem, nodata=np.nan)
    reproject_surface = reproject_surface.reindex()
    dem = dem.reindex()

    ## merge the two rasters
    mask = reproject_surface.isnull() & ~dem.isnull()
    if (value_type == 'elevation'):
        merged_dem = dem.where(mask, reproject_surface)
    elif (value_type == 'depth+') or (value_type == 'depth'):
        
        merged_dem = dem.where(mask, dem - reproject_surface)
    elif (value_type == 'depth-'):
        merged_dem = dem.where(mask, dem + reproject_surface)
    merged_dem.rio.to_raster(merged_file)

    return merged_dem
