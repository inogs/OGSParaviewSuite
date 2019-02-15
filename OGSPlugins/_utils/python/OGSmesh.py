#!/usr/bin/env pvpython
#
# Python class to deal with mesh conversion to ParaView
#
# Arnau Miro, OGS (2018)

import os, argparse
import numpy as np, ctypes as ct

from mpl_toolkits.basemap import Basemap

# Imports from bit.sea python module
from commons.mask import Mask
from commons.submask import SubMask
from basins import V2 as OGS

# Definitions from CTYPES
c_void_p   = ct.c_void_p
c_char     = ct.c_char
c_char_p   = ct.POINTER(ct.c_char)
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
	def __init__(self,maskpath,maskname="meshmask.nc",maptype='merc',lib='libOGS.so'):
		'''
		Class constructor for OGSmesh.

		Inputs:
			> maskpath: Full path to meshmask.
			> maptype:  Kind of map projection (default: merc).
			> lib:      Full path to the OGSmesh.so library.
		'''
		# Class variables
		self.maskpath = maskpath
		self.maskname = maskname
		self.mask     = None     # meshmask will update after being read
		self.map      = maptype

		# Interface with the C functions
		self.OGSlib   = ct.cdll.LoadLibrary(lib)
		
	def readMeshMask(self,fname):
		'''
		Reads the meshmask.nc file and extracts the dimensions of the mesh
		as well as the longitude, latitude and nav_lev vectors.

		Requires NetCDF4 and bit.sea to work.

		Inputs:
			> fname: full path to the meshmask.nc file
		'''
		# Load the meshmask (requires bit.sea)
		# This first mask is for loading the cell centered masks
		self.mask = Mask(fname,zlevelsvar="nav_lev",ylevelsmatvar="gphit", xlevelsmatvar="glamt")
		# This mask is for generating the mesh points
		mask = Mask(fname,zlevelsvar="gdepw",ylevelsmatvar="gphif", xlevelsmatvar="glamf")
		# Set variables
		dims    = list(mask.shape) # We need to add the missing points
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
		return dims, Lon, Lat, nav_lev

	def getMeshResolution(self,dims):
		'''
		Returns the mesh resolution as a string (either "low", "mid" or "high")
		according to the dimensions of the mesh read.

		Inputs:
			> dims: dimensions.
		'''
		if dims == (43, 160, 394):   return "low"
		if dims == (72, 253, 722):   return "mid"
		if dims == (125, 380, 1085): return "high"

		return None

	def generateBasinsMask(self):
		'''
		Generate the basins_mask field where all basins are numbered from 1
		to the number of basins.
		'''
		# Initialize 
		SUBlist     = [ sub.name for sub in OGS.P.basin_list ]
		basins_mask = np.zeros(self.mask.shape)

		# Run for each sub basin
		for sub in SUBlist:
			# Avoid dealing with the whole Mediterranean sea
			if sub == "med": continue
			# Obtain index and basin name
			index = SUBlist.index(sub)
			basin = OGS.P.basin_list[index]
			# Extract the sub mask
			s = SubMask(basin, maskobject=self.mask)
			# Build the basins mask for ParaView
			basins_mask[s.mask] = index + 1

		return basins_mask

	def generateCoastsMask(self):
		'''
		Generate the continetnal shelf mask field where coastal areas 
		are separated from the open sea (at depth 200m).

		FIX AMAL: the intersection with the med mask and mask_at_level
		returns zero for the last value of the mask, even if the depth is less
		than 200. To avoid that, we get the next depth value from 200.
		'''
		dims = self.mask.shape

		# Extract mask at level 200
		# This is all the places that have water at depth = 200 m
		jk_m       = self.mask.getDepthIndex(200.)
		mask200_2D = self.mask.mask[jk_m+1,:,:].copy() # FIX AMAL
		mask200_3D = np.array([mask200_2D for i in xrange(dims[0])])
		
		# Extract mask for mediterranean sea
		# We want all the places that belong to the MED and that are water from 0 to 200 m
		s = SubMask(OGS.P.basin_list[-1], maskobject=self.mask)

		# Define coasts mask
		coasts_mask = np.zeros(dims)
		coasts_mask[~mask200_3D & s.mask] = 1 # Coast
		coasts_mask[ mask200_3D & s.mask] = 2 # Open sea

		return coasts_mask	

	def applyProjection(self,dims,Lon,Lat,Lon0=0.,Lat0=0.):
		'''
		Applies a map projection to the longitude and latitude
		vectors. The kind of map projection is defined in the
		class definition.

		Inputs:
			> dims: dimensions.
			> Lon:  longitude vector.
			> Lat:  latitude vector.
		'''
		# Define map projection
		mproj = Basemap(projection = self.map, \
						lat_0      = Lon0, \
						lon_0      = Lat0, \
                        llcrnrlon  = -5.3, \
                        llcrnrlat  = 28.0, \
                        urcrnrlon  = 37,   \
                        urcrnrlat  = 46.0, \
                        resolution = 'l'
                       )
		# Initialize arrays
		nLon = dims[2]
		nLat = dims[1]
		Lon2Meters = np.zeros((nLon,),np.double)
		Lat2Meters = np.zeros((nLat,),np.double)
		# Perform projection
		for ii in xrange(0,nLon):
			xpt,ypt        = mproj(Lon[60,ii],Lat[60,0]) # FIXED NEW
			Lon2Meters[ii] = xpt
		for jj in xrange(0,nLat):
			xpt,ypt        = mproj(Lon[0,nLon/2],Lat[jj,nLon/2])
			Lat2Meters[jj] = ypt
		# Return
		return Lon2Meters, Lat2Meters

	def OGS(self,fname,Lon2Meters,Lat2Meters,nav_lev,basins_mask,coast_mask):
		'''
		Interface with the C++ OGS class. Create a new object and return
		the pointer to said object.

		Inputs:
			> fname:       Name of the file to write
			> Lon2Meters:  Longitude to meters conversion (npy array)
			> Lat2Meters:  Latitude to meters conversion (npy array)
			> nav_lev:     Depth (npy array)
			> basins_mask: Mask contanining the sub basins
			> coast_mask:  Mask contanining the coasts
		'''
		OGSnew       = self.OGSlib.newOGS
		OGS.argtypes = [c_char_p,c_int,c_int,c_int,c_double_p,c_double_p,c_double_p,c_double_p]
		OGS.restype  = c_void_p

		# Compute sizes of vectors
		nLon = Lon2Meters.shape[0]
		nLat = Lat2Meters.shape[0]
		nLev = nav_lev.shape[0]

		# Return class instannce
		return OGSnew(fname,c_int(nLon),c_int(nLat),c_int(nLev),Lon2Meters.ctypes.data_as(c_double_p),\
					  Lat2Meters.ctypes.data_as(c_double_p), nav_lev.ctypes.data_as(c_double_p),\
					  basins_mask.ctypes.data_as(c_double_p),coast_mask.ctypes.data_as(c_double_p)
					 )

	def OGSwriteMesh(self,OGScls):
		'''
		Wrapper for the C function writeOGSMesh inside the OGSmesh.so library. 
		'''
		writeMesh          = self.OGSlib.OGSWriteMesh
		writeMesh.argtypes = [c_void_p]
		writeMesh.restype  = c_int

		return writeMesh(OGScls)

	def createOGSMesh(self,fname="mesh.ogsmsh"):
		'''
		creates a binary file containing all the mesh information for ParaView as well as
		the basins and coasts masks.

		Inputs:
			> fname: name of the output file to write. It will be written at the same path 
			as the meshmask.nc file (Default: mesh.ogsmsh).
		'''
		# Read the mesh mask
		dims, Lon, Lat, nav_lev = self.readMeshMask(os.path.join(self.maskpath,self.maskname))
		# Obtain the coasts_mask and the basins_mask
		basins_mask = self.generateBasinsMask().ravel()
		coasts_mask = self.generateCoastsMask().ravel()

		# Project latitude and longitude according to map specifics
		Lon2Meters, Lat2Meters = self.applyProjection(dims,Lon,Lat)

		# Create an instance of the OGS class
		OGScls = self.OGS(fname,Lon2Meters,Lat2Meters,nav_lev,basins_mask,coasts_mask);

		# Save into file
		self.OGSwriteMesh(OGScls)

'''
	MAIN

The program runs here without the need to import the class.

'''
if __name__ == '__main__':
	argpar = argparse.ArgumentParser(prog="OGSmesh",description="Convert mesh for ParaView usage.")
	argpar.add_argument('-i','--input',type=str,help='Full path to meshmask directory',required=True,dest='inpath')
	argpar.add_argument('-o','--output',type=str,help='Name of the output mesh file',dest='outfile')
	argpar.add_argument('-m','--map',type=str,help='Projection type (default: merc)',dest='map')

	# parse input arguments
	args=argpar.parse_args()
	if not args.map: args.map = 'merc'

	# Define class instance
	mesh = OGSmesh(args.inpath,args.map)

	# Generate the mesh file
	if not args.outfile:
		mesh.createOGSMesh()
	else:
		mesh.createOGSMesh(args.outfile)