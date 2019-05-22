#!/usr/bin/env python
#
# Conversion of a longitude/latitude pair
# to meters according a desired projection.
#
# Arnau Miro, OGS (2019)

from __future__ import print_function

import numpy as np, argparse
from mpl_toolkits.basemap import Basemap

def lonlat2m(lon,lat,proj='merc',lon0=0.,lat0=0.):
	'''
	Conversion of a longitude/latitude pair to meters 
	according a desired projection using the Basemap
	library.

	This function is provided as a convenience method
	in case the module is imported into ParaView.
	'''
	# Define map projection
	mproj = Basemap(projection = proj, \
					lat_0      = lat0, \
					lon_0      = lon0, \
					llcrnrlon  = -5.3, \
					llcrnrlat  = 28.0, \
					urcrnrlon  = 37,   \
					urcrnrlat  = 46.0, \
					resolution = 'l'
				   )
	# Obtain the coordinates in meters
	xpt,ypt = mproj(lon,lat)
	# Transform the coordinates
	Lon2Meters = xpt if not proj == 'cyl' else 6371e3*np.deg2rad(xpt)
	Lat2Meters = ypt if not proj == 'cyl' else 6371e3*np.deg2rad(ypt)
	# Return
	return Lon2Meters,Lat2Meters	


if __name__ == '__main__':
	argpar = argparse.ArgumentParser(prog="OGSlonlat2m",description="Convert a longitude/latitude pair to meters.")
	argpar.add_argument('--lon', type=float,help='Longitude in degrees', required=True, dest='lon')
	argpar.add_argument('--lat', type=float,help='Latitude in degrees',  required=True, dest='lat')
	argpar.add_argument('-m', '--map', type=str, help='Projection type (default: merc)', dest='map')
	argpar.add_argument('--lon0', type=float,help='Initial longitude in degrees (default: 0)', dest='lon0')
	argpar.add_argument('--lat0', type=float,help='Initial latitude in degrees (default: 0)',  dest='lat0')

	# parse input arguments
	args=argpar.parse_args()
	if not args.map:  args.map  = 'merc'
	if not args.lon0: args.lon0 = 0.
	if not args.lat0: args.lat0 = 0.

	print('x = %f\ny = %f' % lonlat2m(args.lon,args.lat,proj=args.map,lon0=args.lon0,lat0=args.lat0))