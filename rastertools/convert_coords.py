from __future__ import annotations

import numpy as np
from shapely.geometry import Point, LineString
import geopandas as gpd
# import jax.numpy as jnp
# from jax import jit

class CoordConverter(object):
	def __init__(self, cl_geom: LineString):
		self.cl_geom = cl_geom
		self.line_verts = [Point(p) for p in cl_geom.coords]
		self.line_dist = np.array([cl_geom.project(p) for p in self.line_verts])
	
	def xy2sn_coord(self,
					x: float|None = None, 
					y: float|None = None,
					z: float|None = None,
					point: Point|None = None
					) -> tuple[float, float, float|None]:
		"""Covert (x, y, z) to (s, n, z) with a centerline
		The point will be used is a Point object is given.
		Otherwise, (x, y, z) should be specified.

		Parameters
		----------
		cl_geom : LineString
			The geometry of the centerline (flowline).
		x : float | None, optional
			x coordinate of the point, by default None
		y : float | None, optional
			y coordinate of the point, by default None
		z : float | None, optional
			z coordinate of the point, by default None
		point : Point | None, optional
			Point object with (x, y, z), by default None

		Returns
		-------
		tuple[float, float, float|None]
			3D: (s, n, z) or 2D: (s, n, None)
		"""

		if (point == None):
			point = Point(x, y, z) if z else Point(x, y)
		s = self.cl_geom.project(point)
		n = self.cl_geom.distance(point)
		z = point.z if point.has_z else None

		if (s > self.cl_geom.length):
			intp_prev, intp_later = self.line_verts[-2:]
		elif (s <= 0):
			intp_prev, intp_later = self.line_verts[0:2]
		else:
			idx_ls = self.line_dist[(self.line_dist - s) < 0]
			start_idx = (s - idx_ls).argmin()
			intp_prev, intp_later = self.line_verts[start_idx:start_idx+2]

		# # Ring: previous point -> interpolated point -> point
		# lr = LinearRing([intp_prev, intp_later, point])       
		
		# # S-N: cw = right N+, ccw = left N-
		# if lr.is_ccw:
		#     n *= -1

		flow_vec = np.array([intp_later.x - intp_prev.x, intp_later.y - intp_prev.y])
		point_vec = np.array([point.x - intp_prev.x, point.y - intp_prev.y])
		n *= np.sign(np.cross(point_vec, flow_vec)) # sign of z value of cross-product

		return s, n, z

	def xy2sn_point(self,
					x: float|None = None, 
					y: float|None = None,
					z: float|None = None,
					point: Point|None = None
					) -> tuple[float, float, float|None]:
		s, n, z = self.xy2sn_coord(x, y, z, point)
		point = Point(s, n) if z is None else Point(s, n, z)
		return point


	def sn2xy_coord(self,
					s: float|None = None, 
					n: float|None = None,
					z: float|None = None,
					point: Point|None = None
					) -> tuple[float, float, float|None]:
		"""_summary_

		Parameters
		----------
		cl_geom : LineString
			The geometry of the centerline (flowline).
		s : float | None, optional
			s coordinate of the point, by default None
		n : float | None, optional
			n coordinate of the point, by default None
		z : float | None, optional
			z coordinate of the point, by default None
		point : Point | None, optional
			Point object with (s, n, z), by default None

		Returns
		-------
		tuple[float, float, float|None]
			3D: (x, y, z) or 2D: (x, y, None)
		"""
		if point:
			s, n, z = point.coords[0]
		intp = self.cl_geom.interpolate(s)
		if (s > self.cl_geom.length):
			intp_prev = self.line_verts[-1]
		elif (s < 0):
			intp_prev = self.line_verts[0]
		elif (s == 0):
			intp_prev, intp = self.line_verts[0], self.line_verts[1]
		else:
			idx_ls = self.line_dist[(self.line_dist - s) < 0]
			start_idx = (s - idx_ls).argmin()
			intp_prev = self.line_verts[start_idx]
		dist = intp.distance(intp_prev)
		# get unit normal vector of the segment
		vec_N = np.array([intp.y - intp_prev.y, intp_prev.x - intp.x]) / dist
		x, y = intp.coords[0] + vec_N * n
		return x, y, z

	def sn2xy_point(self,
					s: float|None = None, 
					n: float|None = None,
					z: float|None = None,
					point: Point|None = None
					) -> tuple[float, float, float|None]:
		x, y, z = self.sn2xy_coord(s, n, z, point)
		point = Point(x, y) if z is None else Point(x, y, z)
		return point