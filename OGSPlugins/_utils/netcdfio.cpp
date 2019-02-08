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
#include "netcdfio.hpp"

#define MAXVAL 1.e15

namespace NetCDF
{
	/* READNETCDF

		Reads a NetCDF4 file given the name of the file, the name of the variable
		to be read and its size.
	*/
	double *readNetCDF(const char *fname, const char *varname, const int n) {

		int fid, varid, retval;
		double *out;

		// Allocate output variable
		out = new double[n];
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

	/* READNETCDF2F

		Reads a NetCDF4 file given the name of the file, the name of the variable
		to be read and its size.

		Stores the variable in a field structure, therefore it is safe.
	*/
	field::Field<double> readNetCDF2F(const char *fname, const char *varname, const int n) {

		int fid, varid, retval;
		field::Field<double> out = field::Field<double>(n,1);

		// Open file for reading
		if ( nc_open(fname,NC_NOWRITE,&fid) != NC_NOERR ) return out;
		// Get the variable id based on its name
		if ( nc_inq_varid(fid,varname,&varid) != NC_NOERR ) return out;
		// Read the data
		nc_get_var_double(fid,varid,out.data());
		// Close the file
		nc_close(fid);
		// Eliminate the missing variables
		for (int ii=0;ii<n;ii++)
			if(out[ii][0] > MAXVAL) out[ii][0] = 0.;

		// Return
		return out;
	}
}