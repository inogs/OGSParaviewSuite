#!/bin/env python
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


class OGSmesh(object):
	def __init__(self,maskpath,maptype='merc',lib='./libOGSmesh.so'):
		'''
		Class constructor for OGSmesh.

		Inputs:
			> maskpath: Full path to meshmask.
			> maptype:  Kind of map projection (default: merc).
			> lib:      Full path to the OGSmesh.so library.
		'''
		self.maskpath = maskpath
		self.mask     = None     # meshmask will update after being read
		self.map      = maptype
		self.lib      = lib

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
		Generate the coasts_mask field where coastal areas are separated
		from the open sea.

		TODO: deal with annaCoast
		'''
		dims = self.mask.shape
		# Extract mask at level 200
		mask200_2D = self.mask.mask_at_level(200.0)
		mask200_3D = np.zeros(dims,dtype=np.bool)
		for ii in range(dims[0]):
			mask200_3D[ii,:,:] = mask200_2D

		# Extract mask for mediterranean sea
		s = SubMask(OGS.P.basin_list[-1], maskobject=self.mask)

		# Define coasts mask
		coasts_mask = np.zeros(dims)
		coasts_mask[~mask200_3D & s.mask] = 1 # Coast
		coasts_mask[ mask200_3D & s.mask] = 2 # Open sea

		return coasts_mask	

	def applyProjection(self,dims,Lon,Lat):
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
		mproj = Basemap(projection = self.map,
						lat_0      = 0.,   \
						lon_0      = 0.,   \
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

	def writeOGSMesh(self,fname,Lon2Meters,Lat2Meters,nav_lev,basins_mask,coast_mask):
		'''
		Wrapper for the C function writeOGSMesh
		inside the OGSmesh.so library. 

		Inputs:
			> fname:       Name of the file to write
			> Lon2Meters:  Longitude to meters conversion (npy array)
			> Lat2Meters:  Latitude to meters conversion (npy array)
			> nav_lev:     Depth (npy array)
			> basins_mask: Mask contanining the sub basins
			> coast_mask:  Mask contanining the coasts
		'''	
		c_double_p = ct.POINTER(ct.c_double)
		OGSmesh    = ct.cdll.LoadLibrary(self.lib)

		# Compute sizes of vectors
		nLon = np.shape(Lon2Meters)[0]
		nLat = np.shape(Lat2Meters)[0]
		nLev = np.shape(nav_lev)[0]
		
		OGSmesh.writeOGSMesh(fname,ct.c_int(nLon),ct.c_int(nLat),ct.c_int(nLev), \
			Lon2Meters.ctypes.data_as(c_double_p),Lat2Meters.ctypes.data_as(c_double_p), \
			nav_lev.ctypes.data_as(c_double_p),basins_mask.ctypes.data_as(c_double_p), \
			coast_mask.ctypes.data_as(c_double_p))

	def createOGSMesh(self,fname="mesh.ogsmsh"):
		'''
		creates a binary file containing all the mesh information for ParaView as well as
		the basins and coasts masks.

		Inputs:
			> fname: name of the output file to write. It will be written at the same path 
			as the meshmask.nc file (Default: mesh.ogsmsh).
		'''
		# Read the mesh mask
		dims, Lon, Lat, nav_lev = self.readMeshMask(os.path.join(self.maskpath,"meshmask.nc"))
		# Obtain the coasts_mask and the basins_mask
		basins_mask = self.generateBasinsMask().ravel()
		coasts_mask = self.generateCoastsMask().ravel()

		# Project latitude and longitude according to map specifics
		Lon2Meters, Lat2Meters = self.applyProjection(dims,Lon,Lat)

		# Save into file
		self.writeOGSMesh(os.path.join(self.maskpath,fname),
						  Lon2Meters,Lat2Meters,nav_lev.astype(np.double),
			              basins_mask.astype(np.double),coasts_mask.astype(np.double)
			             )

'''
	MAIN

The program runs here without the need to import the class.

'''
if __name__ == '__main__':
	argpar = argparse.ArgumentParser(prog="OGSmesh",description="Convert mesh for ParaView usage.")
	argpar.add_argument('-i','--input',type=str,help='Full path to meshmask directory',required=True,dest='inpath')
	argpar.add_argument('-o','--output',type=str,help='Name of the output mesh file',dest='outfile')
	argpar.add_argument('-m','--map',type=str,help='Projection type (default: merc)',dest='map')
	argpar.add_argument('-l','--lib',type=str,help='Path to the libOGSmesh.so library',dest='lib')

	# parse input arguments
	args=argpar.parse_args()
	if not args.map: args.map = 'merc'
	if not args.lib: args.lib = './libOGSmesh.so'

	# Define class instance
	mesh = OGSmesh(args.inpath,args.map,args.lib)

	# Generate the mesh file
	if not args.outfile:
		mesh.createOGSMesh()
	else:
		mesh.createOGSMesh(args.outfile)