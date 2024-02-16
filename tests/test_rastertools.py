import pytest
from rastertools import *
import geopandas as gpd
import numpy as np
import rasterio as rio
import rioxarray
from rasterio.transform import Affine
import os

# def test_interpolate_points_to_raster(mock_mesh_file, mock_surface_file):


#     # Call the interpolate_points_to_raster function with the mock files
#     output_file = mock_surface_file
#     interpolate_points_to_raster(mock_mesh_file, surface_file=output_file, resolution=1)

#     # Open the output file and check the interpolated values
#     interpolated = rioxarray.open_rasterio(mock_mesh_file, masked=True)
#     expected_values = np.array([[0, 1, 2], [1, 2, 3], [2, 3, 4]], dtype=np.float32)
#     np.testing.assert_array_equal(interpolated.data, expected_values)

#     # Clean up the output file
#     os.remove(output_file)


