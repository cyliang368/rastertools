#%%
import urllib, requests
import json
import geopandas as gpd
import rasterio as rio
from rasterio.mask import mask as raster_mask
import matplotlib.pyplot as plt
import os
import pyproj
from shapely.geometry import *
from osgeo import gdal, osr

class RasterTools(object):
	def __init__(self, root_dir=None, centerline_file=None, DEM_file=None, DEM_path=None):
		self.root_dir = root_dir if root_dir else './temp'
		if (not os.path.exists(self.root_dir)):
			os.makedirs(self.root_dir)
		
		if (isinstance(centerline_file, str)):
			self.bd_gdf = gpd.read_file(centerline_file)
			self.crs_org = self.bd_gdf.crs
		elif (isinstance(centerline_file, gpd.GeoDataFrame)):
			self.bd_gdf = centerline_file
			self.crs_org = self.bd_gdf.crs
		else:
			self.bd_gdf = None
			self.crs_org = None
		
		self.DEM_file = DEM_file
		# if (self.DEM_path):
		# 	self.DEM = rio.open(DEM_path)
		
	def getDEMbyPolygon(self, centerline_file=None, output_dir = None, output_file= None, format ='GeoTIFF', grid_size=r'1/3%20arcsecond'):
		if (centerline_file is None):
			poly_gdf = self.bd_gdf
		elif (isinstance(centerline_file, str)):
			poly_gdf = gpd.read_file(centerline_file)
		elif (isinstance(centerline_file, gpd.GeoDataFrame)):
			poly_gdf = centerline_file
		poly_gdf = poly_gdf.to_crs(4326)
	
		## get boundary box from the polygon geometry in the file
		bounds = poly_gdf.loc[0, 'geometry'].bounds
		bounds_str = str(bounds)[1:-1].replace(' ', '')
		bounds_str = bounds_str

		url = r'https://tnmaccess.nationalmap.gov/api/v1/products?datasets=National%20Elevation%20Dataset%20(NED)%20'
		url += r'1/3%20arc-second'
		url += r'&&bbox=' + bounds_str
		url += r'&prodFormats=' + format
		url += r'&dateType=' + 'Publication'
		url += r'&start=' + '2006-01-01'
		url += r'&end=' + '2022-01-01'

		self.query = requests.get(url, verify=False) 
		self.query_json_dic = self.query.json()

		downloadURL = self.query_json_dic['items'][0]['downloadURL']
		if (output_file):
			output_path = os.path.join(output_dir, output_file)
		else:
			output_dir = output_dir if output_dir else self.root_dir
			output_path = os.path.join(output_dir, os.path.basename(downloadURL))

		output_path, _ = urllib.request.urlretrieve(downloadURL, output_path)
		self.DEM_file = output_path

	def clipDEMbyPolygonExtent(self, polygon_file=None, input_file=None, output_file=None, buffer_ratio=1.0, buffer_pixels=0):
		if (polygon_file is None):
			poly_gdf = self.bd_gdf
		elif (isinstance(polygon_file, str)):
			poly_gdf = gpd.read_file(polygon_file)
		elif (isinstance(polygon_file, gpd.GeoDataFrame)):
			poly_gdf = polygon_file
		
		## get source DEM resoluations
		input_file = input_file if input_file else self.DEM_file
		ds = gdal.Open(input_file)
		gtf = ds.GetGeoTransform()
		res_x, res_y = abs(gtf[1]), abs(gtf[5])

		## get source DEM crs and project polygon to this crs
		prj = ds.GetProjection()
		srs = osr.SpatialReference(wkt=prj)
		crs = srs.GetAttrValue("AUTHORITY", 1)
		poly_gdf.to_crs(epsg=crs, inplace=True)

		## get the boundary box for clip
		xmin, ymin, xmax, ymax = poly_gdf.loc[0, 'geometry'].bounds
		xmax += res_x
		ymin -= res_y
		if (buffer_pixels > 0):
			xmin -= res_x * buffer_pixels / 2
			xmax += res_x * buffer_pixels / 2
			ymin -= res_y * buffer_pixels / 2
			ymax += res_y * buffer_pixels / 2

		## clip the raster and save as the outputfile
		ds = gdal.Translate(output_file, ds, projWin=[xmin, ymax, xmax, ymin])
		ds = None

	def reprojectRaster(self, input_file, output_file, crs=None):
		if (crs):
			dest_crs = crs
		else:
			dest_crs = self.crs_org
		
		input_raster = gdal.Open(input_file)
		warp = gdal.Warp(output_file, input_raster, dstSRS=dest_crs)
		warp = None
		
	
if __name__ == '__main__':
	clipped_file = './temp/USGS_13_n30w096_20211103_clipped.tif'
	reprojected_file = './temp/USGS_13_n30w096_20211103_clipped_proj.tif'

	bd_file = '../XS_processing/xs_tools/Brazos/input/boundary.shp'
	# rt = RasterTools(centerline_file=bd_file, DEM_file='./temp/USGS_13_n30w096_20211103.tif')
	rt = RasterTools(centerline_file=bd_file)
	rt.getDEMbyPolygon()
	rt.clipDEMbyPolygonExtent(output_file=clipped_file, buffer_pixels=100)
	rt.reprojectRaster(input_file=clipped_file, output_file=reprojected_file)

# %%
