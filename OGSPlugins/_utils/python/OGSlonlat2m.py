#!/usr/bin/env python
#
# Conversion of a longitude/latitude pair
# to meters according a desired projection.
#
# (c) OGS, Arnau Miro (2019)
from __future__ import print_function

import numpy as np, argparse
from OGSmesh import OGSmesh

if __name__ == '__main__':
	argpar = argparse.ArgumentParser(prog="OGSlonlat2m",description="Convert a longitude/latitude pair to meters.")
	argpar.add_argument('--lon',type=float,help='Longitude in degrees', required=True, dest='lon')
	argpar.add_argument('--lat',type=float,help='Latitude in degrees',  required=True, dest='lat')
	argpar.add_argument('-p','--proj',type=str,help='Projection type (default: Mercator)',dest='proj')
	argpar.add_argument('--with-kwags',action='store_true',help='Force launch with predefined kwargs',dest='withkwargs')

	# parse input arguments
	args=argpar.parse_args()
	if not args.proj:  args.proj  = 'Mercator'

	mesh = OGSmesh('./')
	projkwargs = {'false_easting':989634.3811336625,'false_northing':-3512473.95569} if args.withkwargs else dict()

	lon2m,lat2m = mesh.applyProjection(args.proj,np.array([args.lon]),np.array([args.lat]),**projkwargs)

	print( 'x = %f\ny = %f' % (lon2m,lat2m) )