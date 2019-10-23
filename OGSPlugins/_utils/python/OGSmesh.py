#!/usr/bin/env pvpython
#
# Python class to deal with mesh conversion to ParaView
#
# (c) OGS, Arnau Miro (2018)
from __future__ import print_function

import os, sys, argparse
import numpy as np, ctypes as ct
import cartopy.crs as ccrs, cartopy.feature as cfeat, shapely.geometry as sgeom

# Imports from bit.sea python module
from commons.mask import Mask
from commons.submask import SubMask
from basins import V2 as OGS

# Definitions from CTYPES
c_void_p   = ct.c_void_p
c_char     = ct.c_char
c_char_p   = ct.POINTER(ct.c_char)
c_uint8    = ct.c_uint8
c_uint8_p  = ct.POINTER(ct.c_uint8)
c_int      = ct.c_int
c_int_p    = ct.POINTER(ct.c_int)
c_double   = ct.c_double
c_double_p = ct.POINTER(ct.c_double)

class OGSmesh(object):
	'''
	Python wrapper to the OGS C++ class that lets us interface with the
	OGS mesh (.ogsmsh) file.

	This class needs the libOGS.so that is generally deployed along with
	this class by the deployment scripts of the OGSParaView Suite.
	'''
	def __init__(self,maskpath,maskname="meshmask.nc",res='simple',lib='libOGS'):
		'''
		Class constructor for OGSmesh.

		Inputs:
			> maskpath   : Full path to the mask file.
			> maskname   : Name of the mask file (default: meshmask.nc)
			> projection : Kind of map projection (default: merc).
			> res        : Resolution for basemap (default: l)
			> lib        : Full path to the OGSmesh.so library.
		'''
		# Class variables
		self._maskpath    = maskpath
		self._maskname    = maskname
		self._res         = res

		# Class properties
		self._mask        = None
		self._basins_mask = []
		self._coasts_mask = []
		self._land_mask   = []

		# Interface with the C functions
		lib += '.dylib' if sys.platform == 'darwin' else '.so'
		self._OGSlib   = ct.cdll.LoadLibrary(lib)

	@property
	def mask(self):
		if self._mask == None:
			_, _, _, _, self._mask = self.readMeshMask(os.path.join(self._maskpath,self._maskname))
		return self._mask
	@mask.setter
	def mask(self,maskfile):
		_, _, _, _, self._mask = self.readMeshMask(maskfile)

	@property
	def basins_mask(self):
		if self._basins_mask == []:
			self._basins_mask = self.generateBasinsMask(self.mask)
		return self._basins_mask
	@basins_mask.setter
	def basins_mask(self,basins_mask):
		self._basins_mask = basins_mask

	@property
	def coasts_mask(self):
		if self._coasts_mask == []:
			self._coasts_mask = self.generateCoastsMask(self.mask)
		return self._coasts_mask
	@coasts_mask.setter
	def coasts_mask(self,coasts_mask):
		self._coasts_mask = coasts_mask

	@property
	def land_mask(self):
		if self._land_mask == []:
			self._land_mask = self.generateLandMask(self._res,self.mask,self.coasts_mask,self.is_land)
		return self._land_mask
	@land_mask.setter
	def land_mask(self,land_mask):
		self._land_mask = land_mask
		
	@staticmethod	
	def readMeshMask(fname):
		'''
		Reads the meshmask.nc file and extracts the dimensions of the mesh
		as well as the longitude, latitude and nav_lev vectors.

		Requires NetCDF4 and bit.sea to work.

		Inputs:
			> fname: full path to the meshmask.nc file
		'''
		# Load the meshmask (requires bit.sea)
		# This first mask is for loading the cell centered masks
		mask1 = Mask(fname,zlevelsvar="nav_lev",ylevelsmatvar="gphit", xlevelsmatvar="glamt")
		# This mask is for generating the mesh points
		mask = Mask(fname,zlevelsvar="gdepw",ylevelsmatvar="gphif", xlevelsmatvar="glamf")
		# Set variables
		dims     = list(mask.shape) # We need to add the missing points
		dims[0] += 1
		dims[1] += 1
		dims[2] += 1
		# Longitudinal coordinates, add one in the right
		Lon = np.insert(mask.xlevels,0,mask.xlevels[0,:] - (mask.xlevels[0,:] - mask.xlevels[1,:])/2.,axis=0)
		Lon = np.insert(Lon,0,Lon[:,0] - (Lon[:,0] - Lon[:,1])/2.,axis=1)
		# Latitudinal coordinates, add one in the right
		Lat = np.insert(mask.ylevels,0,mask.ylevels[0,:] + (mask.ylevels[0,:] - mask.ylevels[1,:])/2.,axis=0)
		Lat = np.insert(Lat,0,Lat[:,0] - (Lat[:,0] - Lat[:,1])/2.,axis=1)
		# Depth coordinates, add the one in the left
		nav_lev = np.append(mask.zlevels,mask.zlevels[-1] + (mask.zlevels[-1] - mask.zlevels[-2])/2.)

		# Return
		return dims, Lon, Lat, nav_lev, mask1

	@staticmethod
	def generateBasinsMask(mask):
		'''
		Generate the basins_mask field where all basins are numbered from 1
		to the number of basins.

		Inputs:
			> mask : the mask object
		'''
		return np.array([SubMask(sub, maskobject=mask).mask.ravel() for sub in OGS.P.basin_list[:-1]],dtype=c_uint8).T

	@staticmethod
	def generateCoastsMask(mask):
		'''
		Generate the continental shelf mask field where coastal areas 
		are separated from the open sea (at depth 200m).

		FIX AMAL: the intersection with the med mask and mask_at_level
		returns zero for the last value of the mask, even if the depth is less
		than 200. To avoid that, we get the next depth value from 200.

		Inputs:
			> mask : the mask object
		'''
		dims = mask.shape

		# Extract mask at level 200
		# This is all the places that have water at depth = 200 m
		jk_m       = mask.getDepthIndex(200.)
		mask200_2D = mask.mask[jk_m+1,:,:].copy() # FIX AMAL
		mask200_3D = np.array([mask200_2D for i in range(dims[0])])
		
		# Extract mask for mediterranean sea
		# We want all the places that belong to the MED and that are water from 0 to 200 m
		s = SubMask(OGS.P.basin_list[-1], maskobject=mask)

		# Define coasts mask
		coasts_mask = np.zeros(dims,dtype=c_uint8)
		coasts_mask[~mask200_3D & s.mask] = 1 # Coast
		coasts_mask[ mask200_3D & s.mask] = 2 # Open sea

		return coasts_mask

	@staticmethod
	def generateLandMask(res,mask,coasts_mask,is_land_fun):
		'''
		Generate the mask field that shows true when a point is in
		a land area and false when it is on the water. This mask accounts
		for the atlantic buffer, the gulf of Bizcaya and the black sea so
		that the continents are easily recognizable.

		Inputs:
			> res         : the resolution for the land ('10m', '50m', '110m' or 'simple')
			> mask        : a mask object
			> coasts_mask : the coasts mask
			> is_land_fun : function that returns if a point is land or not
		'''
		# Create a containing array
		dims = mask.shape
		land_mask = np.zeros(dims,dtype=c_uint8)

		# Anything that is not water will be land
		land_mask[coasts_mask == 0] = 1

		if not res == 'simple':
			# Load the land feature
			g_land = [g for g in cfeat.NaturalEarthFeature('physical', 'land', res).geometries()]
			ig = 0

			# Loop the surface level and decide whether the point
			# belongs to land or water.
			for jj in range(dims[1]):
				for ii in range(dims[2]):
					# Abort the points that are in the water to capture the battimetry
					if not land_mask[0,jj,ii]: continue

					# The points of the meshmask are in degrees, which is what cartopy needs
					xpt = mask.xlevels[jj,ii]
					ypt = mask.ylevels[jj,ii]

					# Update the land mask
					island,ig = is_land_fun(xpt,ypt,g_land,ig)
					land_mask[:,jj,ii] = 1 if island else 0

		return land_mask

	@staticmethod
	def applyProjection(proj,Lon,Lat,**kwargs):
		'''
		Applies a map projection to the longitude and latitude
		vectors. The kind of map projection is defined in the
		class definition.

		Inputs:
			> proj    : kind of map projection.
			> Lon     : longitude vector.
			> Lat     : latitude vector.
			> **kargs : arguments when creating the projection
		'''			
		# Load the PlateCarree which lets us transform degrees to meters
		crs_degs = ccrs.PlateCarree()
		crs_proj = None

		# Recover the kind of projection from cartopy
		try:
			# Projection can either be a name or a epsg nomenclature
			# check if we are dealing with epsg
			if (proj.find('epsg') >= 0):
				# We have an epsg code in the format epsg:<code>
				# pyepsg needs to be installed for that to work
				crs_proj = ccrs.epsg(code=int(proj.split('_')[1]))
			else:
				# We have a normal projection name
				# Let's first check for two particular names
				if (proj.lower() == 'cylindrical'):
					# Use epsg4087 : WGS 84 / World Equidistant Cylindrical (https://epsg.io/4087)
					crs_proj = ccrs.epsg(4087)
				elif (proj.lower() == 'google'):
					# Google Mercator projection (https://scitools.org.uk/cartopy/docs/latest/cartopy_outline.html?highlight=google_mercator)
					crs_proj = ccrs.GOOGLE_MERCATOR
				else:
					# Use the name to generate a projection
					crs_proj = getattr(ccrs,proj)(**kwargs)
		except:
			raise ValueError('Problems with projection %s!' % proj)

		# At this point we should have a projection and 
		# we should be ready to perform the projection

		if (Lon.shape[0] > 1):
			out = crs_proj.transform_points(crs_degs,Lon,Lat)

			# Obtain Lon2Meters and Lat2Meters to generate the rectilinear grid.
			Lon2Meters = np.sort(out[60,:,0])              # Longitude 2 meters
			Lat2Meters = np.sort(out[:,out.shape[1]//2,1]) # Latitude  2 meters
		else:
			out = crs_proj.transform_point(Lon,Lat,crs_degs)
			Lon2Meters = out[0]
			Lat2Meters = out[1]

		return Lon2Meters, Lat2Meters

	@staticmethod
	def is_land(lon,lat,geoms,ig=0):
		'''
		Given  a point in degrees return if it belongs to land

		Inputs:
			> lon   : longitude of the point (in degrees)
			> lat   : latitude of the point (in degrees)
			> geoms : geometries of the land areas
			> ig    : previous iteration point, for a faster algorithm
		'''
		point = sgeom.Point(lon,lat)

		# check if the point is contained on the same ig element
		if geoms[ig].contains(point):
			return True, ig

		# scan from ig to the end
		lg = len(geoms)
		for ii in range(ig,lg):
			if geoms[ii].contains(point):
				return True, ii

		# scan from 0 to ig
		for ii in range(0,ig):
			if geoms[ii].contains(point):
				return True, ii

		# the point is not land
		return False, ig

	def OGSwriteMesh(self,fname,wrkdir,Lon2Meters,Lat2Meters,nav_lev,basins_mask,coast_mask,land_mask):
		'''
		Wrapper for the C function writeOGSMesh inside the OGSmesh.so library. 
		'''
		writeMesh          = self._OGSlib.OGSWriteMesh
		writeMesh.argtypes = [c_char_p,c_char_p,c_int,c_int,c_int,c_double_p,c_double_p,c_double_p,c_uint8_p,c_uint8_p,c_uint8_p]
		writeMesh.restype  = c_int

		# Compute sizes of vectors
		nLon = Lon2Meters.shape[0]
		nLat = Lat2Meters.shape[0]
		nLev = nav_lev.shape[0]

		# Return class instance
		return writeMesh(fname.encode('utf-8'),wrkdir.encode('utf-8'),c_int(nLon),c_int(nLat),c_int(nLev),\
						 Lon2Meters.ctypes.data_as(c_double_p),\
					     Lat2Meters.ctypes.data_as(c_double_p),\
					     nav_lev.ctypes.data_as(c_double_p),\
					     basins_mask.ctypes.data_as(c_uint8_p),\
					     coast_mask.ctypes.data_as(c_uint8_p),\
					     land_mask.ctypes.data_as(c_uint8_p)
					    )

	def createOGSMesh(self,fname="mesh.ogsmsh",path=".",proj='Mercator',projkwargs=dict()):
		'''
		creates a binary file containing all the mesh information for ParaView as well as
		the basins and coasts masks.

		Inputs:
			> fname: name of the output file to write. It will be written at the same path 
			as the meshmask.nc file (Default: mesh.ogsmsh).
		'''
		# Read the mesh mask
		dims, Lon, Lat, nav_lev, self._mask = self.readMeshMask(os.path.join(self._maskpath,self._maskname))

		# Project latitude and longitude according to map specifics
		Lon2Meters, Lat2Meters = self.applyProjection(proj,Lon,Lat,**projkwargs)

		# Save into file
		self.OGSwriteMesh(fname,path,Lon2Meters,Lat2Meters,nav_lev,
			self.basins_mask.ravel(),self.coasts_mask.ravel(),self.land_mask.ravel());

'''
	MAIN

The program runs here without the need to import the class.

'''
if __name__ == '__main__':
	argpar = argparse.ArgumentParser(prog="OGSmesh",description="Convert mesh for ParaView usage.")
	argpar.add_argument('-i','--input',type=str,help='Full path to meshmask directory',required=True,dest='inpath')
	argpar.add_argument('-o','--output',type=str,help='Name of the output mesh file',dest='outfile')
	argpar.add_argument('-p','--proj',type=str,help='Projection type (default: Mercator)',dest='proj')
	argpar.add_argument('-r','--res',type=str,help='Map resolution (default: simple)',dest='res')

	# parse input arguments
	args=argpar.parse_args()
	if not args.outfile: args.outfile = 'mesh.ogsmsh'
	if not args.proj:    args.proj    = 'Mercator'
	if not args.res:     args.res     = 'simple'

	# Define class instance
	mesh = OGSmesh(args.inpath,res=args.res)
	projkwargs = {'false_easting':989634.3811336625,'false_northing':-3512473.95569}

	# Generate the mesh file
	mesh.createOGSMesh(fname=args.outfile,proj=args.proj,projkwargs=projkwargs)