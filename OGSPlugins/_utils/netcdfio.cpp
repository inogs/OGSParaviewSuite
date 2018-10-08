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

double *readNetCDF(const char *fname, const char *varname, int n) {

	int fid, varid, retval;
	double *out;

	// Allocate output variable
	out = (double*)malloc(n*sizeof(double));
	// Open file for reading
	NetCDF::nc_open(fname,NC_NOWRITE,&fid);
	// Get the variable id based on its name
	NetCDF::nc_inq_varid(fid,varname,&varid);
	// Read the data
	NetCDF::nc_get_var_double(fid,varid,out);
	// Close the file
	NetCDF::nc_close(fid);
	// Eliminate the missing variables
	for (int ii=0;ii<n;ii++)
		if(out[ii] > MAXVAL) out[ii] = 0.;

	// Return
	return out;
}