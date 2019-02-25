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

		int fid, varid;
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
	int readNetCDF2F(const char *fname, const char *varname, field::Field<double> &f) {

		int fid, varid;

		// Open file for reading
		if ( nc_open(fname,NC_NOWRITE,&fid) != NC_NOERR )   return NETCDF_ERR;
		// Get the variable id based on its name
		if ( nc_inq_varid(fid,varname,&varid) != NC_NOERR ) return NETCDF_ERR;
		
		// Read the data
		nc_get_var_double(fid,varid,f.data());
		// Close the file
		nc_close(fid);
		
		// Eliminate the missing variables
		field::Field<double>::iterator iter;
		for (iter = f.begin(); iter != f.end(); ++iter)
			if (iter[0] > MAXVAL) iter[0] = 0.;

		return NETCDF_OK;
	}

	int readNetCDF2F(const char *fname, const char *varname, field::Field<float> &f) {

		int fid, varid;

		// Open file for reading
		if ( nc_open(fname,NC_NOWRITE,&fid) != NC_NOERR )   return NETCDF_ERR;
		// Get the variable id based on its name
		if ( nc_inq_varid(fid,varname,&varid) != NC_NOERR ) return NETCDF_ERR;
		
		// Read the data
		nc_get_var_float(fid,varid,f.data());
		// Close the file
		nc_close(fid);
		
		// Eliminate the missing variables
		field::Field<float>::iterator iter;
		for (iter = f.begin(); iter != f.end(); ++iter)
			if (iter[0] > MAXVAL) iter[0] = 0.;

		return NETCDF_OK;
	}

	/* READNETCDF2F3

		Reads a NetCDF4 file given the name of the file, the name of the variable
		to be read and its size (vector field).

		Stores the variable in a field structure, therefore it is safe.
	*/
	int readNetCDF2F3(const char *fname, const char *vname1, const char *vname2, 
		const char *vname3, field::Field<double> &f) {

		int fid, varid1, varid2, varid3, n = f.get_n();
		double *u, *v, *w;

		// Allocate
		u = new double[n]; v = new double[n]; w = new double[n];

		// Open file for reading
		if ( nc_open(fname,NC_NOWRITE,&fid) != NC_NOERR )   return NETCDF_ERR;
		// Get the variable id based on its name
		if ( nc_inq_varid(fid,vname1,&varid1) != NC_NOERR ) return NETCDF_ERR;
		if ( nc_inq_varid(fid,vname2,&varid2) != NC_NOERR ) return NETCDF_ERR;
		if ( nc_inq_varid(fid,vname3,&varid3) != NC_NOERR ) return NETCDF_ERR;
		
		// Read the data
		nc_get_var_double(fid,varid1,u);
		nc_get_var_double(fid,varid2,v);
		nc_get_var_double(fid,varid3,w);
		// Close the file
		nc_close(fid);

		// Set field and eliminate the missing variables
		field::Field<double>::iterator iter;
		for (iter = f.begin(); iter != f.end(); ++iter) {
			iter[0] = u[iter.ind()]; if(iter[0] > MAXVAL) iter[0] = 0.;
			iter[1] = v[iter.ind()]; if(iter[1] > MAXVAL) iter[1] = 0.;
			iter[2] = w[iter.ind()]; if(iter[2] > MAXVAL) iter[2] = 0.;
		}

		// Return
		delete [] u; delete [] v; delete [] w;
		return NETCDF_OK;
	}
	int readNetCDF2F3(const char *fname, const char *vname1, const char *vname2, 
		const char *vname3, field::Field<float> &f) {

		int fid, varid1, varid2, varid3, n = f.get_n();
		float *u, *v, *w;

		// Allocate
		u = new float[n]; v = new float[n]; w = new float[n];

		// Open file for reading
		if ( nc_open(fname,NC_NOWRITE,&fid) != NC_NOERR )   return NETCDF_ERR;
		// Get the variable id based on its name
		if ( nc_inq_varid(fid,vname1,&varid1) != NC_NOERR ) return NETCDF_ERR;
		if ( nc_inq_varid(fid,vname2,&varid2) != NC_NOERR ) return NETCDF_ERR;
		if ( nc_inq_varid(fid,vname3,&varid3) != NC_NOERR ) return NETCDF_ERR;
		
		// Read the data
		nc_get_var_float(fid,varid1,u);
		nc_get_var_float(fid,varid2,v);
		nc_get_var_float(fid,varid3,w);
		// Close the file
		nc_close(fid);

		// Set field and eliminate the missing variables
		field::Field<float>::iterator iter;
		for (iter = f.begin(); iter != f.end(); ++iter) {
			iter[0] = u[iter.ind()]; if(iter[0] > MAXVAL) iter[0] = 0.;
			iter[1] = v[iter.ind()]; if(iter[1] > MAXVAL) iter[1] = 0.;
			iter[2] = w[iter.ind()]; if(iter[2] > MAXVAL) iter[2] = 0.;
		}

		// Return
		delete [] u; delete [] v; delete [] w;
		return NETCDF_OK;
	}
}