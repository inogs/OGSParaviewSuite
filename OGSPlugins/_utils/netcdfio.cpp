/*=========================================================================

  Program:   Utilities
  Module:    netcdfio.cpp

  This module handles the reading and eventual writing of NetCDF files
  for the OGS Paraview Suite

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

// Include NetCDF functions
// The ParaView NetCDF library is linked here instead of an external one
// so as not to incur in compilation duplicates
#include "vtknetcdf/include/netcdf.h"

#define MAXVAL 1.e15

/* READNETCDF

	Reads a NetCDF4 file given the name of the file, the name of the variable
	to be read and its size.
*/
double *readNetCDF(const char *fname, const char *varname, int n) {

	int fid, varid, retval;
	double *out;

	// Allocate output variable
	out = (double*)malloc(n*sizeof(double));
	// Open file for reading
	if ( nc_open(fname,NC_NOWRITE,&fid) != NC_NOERR )
		return NULL;
	// Get the variable id based on its name
	if ( nc_inq_varid(fid,varname,&varid) != NC_NOERR )
		return NULL;
	// Read the data
	nc_get_var_double(fid,varid,out);
	// Close the file
	nc_close(fid);
	// Eliminate the missing variables
	for (int ii=0;ii<n;ii++)
		if(out[ii] > MAXVAL) out[ii] = 0.;

	// Return
	return out;
}