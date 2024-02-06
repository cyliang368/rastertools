# rastertools

rastertools is a Python package that provides functionality for processing rasters for riverine analyses.


<!-- ## Installation

You can install cstools using ``pip``:

    $ pip install cstools

or using ``conda`` (``mamba``):

    $ conda install -c conda-forge cstools -->


## Quick start

The following example shows how to generate raster from points and intergrate mesh with DEM.

``interpolate_points_to_raster`` in ``rastertools`` creates a raster from very dense points (mesh). A cubic polynomial interpolation is used here. Note that ``boundary_file`` is necessary if we plan to integrate the surface resulting from interpolation with a DEM. The output raster from the algorithm covers the convex hull of input points. Therefore, a boundary may be needed to cut out the unnecessary areas.  ``resolution`` here is used to specify the size of grid if change is needed. The unit of the ``resolution`` is meter, so every layer will be reprojected to a UTM Zone CRS using SI unit.

``integrate_surface_with_dem`` in ``rastertools`` merges ``surface_file`` and ``dem_file`` together. The ``merged_file`` can be considered a output path. 

The meaning of z values in the surface should be specified by ``value_type``: ``elevation``, ``depth+``, ``depth``, and ``depth-``.
``elevation``: the elevation in the given datum.
``depth+`` and ``depth``: the depth from the water surface, the point below water has a positive value
``depth+`` and ``depth``: the depth from the water surface, the point below water has a negative value


```python
    from rastertools import *

    folder = "../data/brazos"
    input_folder = folder + "/input"
    output_folder = folder + "/output"

    mesh_file = input_folder + "/mesh_points.geojson"
    boundary_file = input_folder + "/boundary.geojson"
    dem_file = input_folder + "/dem_clip.tif"
    surface_file = output_folder + "/surface.tif"
    merged_file = output_folder + "/merged.tif"
    
    value_type = 'depth+'

    interpolate_points_to_raster(mesh_file=mesh_file, 
                                surface_file=surface_file, 
                                boundary_file=boundary_file, 
                                resolution=5)
    
    
    integrate_surface_with_dem(surface_file=surface_file,
                            dem_file=dem_file,
                            merged_file=merged_file, 
                            value_type=value_type)
```
