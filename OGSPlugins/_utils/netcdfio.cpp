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

#include <algorithm>

#ifdef __linux__
// Include OpenMP when working with GCC
#include <omp.h>
#define OMP_NUM_THREADS omp_get_num_threads()
#define OMP_THREAD_NUM  omp_get_thread_num()
#else
#define OMP_NUM_THREADS 1
#define OMP_THREAD_NUM  0
#endif

#define PNTIND(ii,jj,kk,nx,ny) ( (nx)*(ny)*(kk) + (nx)*(jj) + (ii) )

namespace NetCDF
{
	/* NETCDF API

		Helper functions for opening, creating, reading and 
		writing NetCDF files.
	*/
	int NetCDFCreate(const char *fname, int &fid) {
		// Create the file. The NC_CLOBBER parameter tells netCDF to overwrite 
		// this file, if it already exists.
		if ( nc_create(fname,NC_CLOBBER,&fid) != NC_NOERR ) return NETCDF_ERR;
		if ( nc_put_att_text(fid, NC_GLOBAL, "Convenctions", strlen("COARDS"), "COARDS") != NC_NOERR ) return NETCDF_ERR;

		return NETCDF_OK;
	}
	int NetCDFDefXY(int fid, int &id, int &vid, const char *dimname, const char *varname, nc_type xtype, int n,
		const char *long_name, const char *units) {
		// Define dimension and variable
		if ( nc_def_dim(fid,dimname,n,&id)            != NC_NOERR ) return NETCDF_ERR;
		if ( nc_def_var(fid,varname,xtype,1,&id,&vid) != NC_NOERR ) return NETCDF_ERR;
		// Define units attributes for the variables
		if ( nc_put_att_text(fid,id,"units",strlen(units),units)             != NC_NOERR ) return NETCDF_ERR;
		if ( nc_put_att_text(fid,id,"long_name",strlen(long_name),long_name) != NC_NOERR ) return NETCDF_ERR;

		return NETCDF_OK;
	}
	int NetCDFDefZ(int fid, int &id, int &vid, const char *dimname, const char *varname, nc_type xtype, int n,
		const char *units, const char *positive) {
		// Define dimension and variable
		if ( nc_def_dim(fid,dimname,n,&id)            != NC_NOERR ) return NETCDF_ERR;
		if ( nc_def_var(fid,varname,xtype,1,&id,&vid) != NC_NOERR ) return NETCDF_ERR;
		// Define units attributes for the variables
		if ( nc_put_att_text(fid,id,"units",strlen(units),units)          != NC_NOERR ) return NETCDF_ERR;
		if ( nc_put_att_text(fid,id,"positive",strlen(positive),positive) != NC_NOERR ) return NETCDF_ERR;

		return NETCDF_OK;
	}
	int NetCDFDefT(int fid, int &id, const char *dimname) {
		if ( nc_def_dim(fid,dimname,NC_UNLIMITED,&id) != NC_NOERR ) return NETCDF_ERR;

		return NETCDF_OK;
	}
	int NetCDFDefVar(int fid, int &id, const char *varname, nc_type xtype, int ndim, int dims[]) {
		if ( nc_def_var(fid,varname,xtype,ndim,dims,&id) != NC_NOERR ) return NETCDF_ERR;
		if ( nc_put_att_text(fid,id,"missing_value",strlen("1.e+20f"),"1.e+20f") != NC_NOERR ) return NETCDF_ERR;

		return NETCDF_OK;
	}
	int NetCDFPutDim(int fid, int id, double *data) {
		if ( nc_put_var_double(fid,id,data) != NC_NOERR ) return NETCDF_ERR;

		return NETCDF_OK;
	}
	int NetCDFPutDim(int fid, int id, float *data) {
		if ( nc_put_var_float(fid,id,data) != NC_NOERR ) return NETCDF_ERR;

		return NETCDF_OK;
	}
	int NetCDFPutVar(int fid, int id, int dims, double *data) {
		// These settings tell netcdf to write one timestep of data. (The
		// setting of start[0] inside the loop below tells netCDF which
		// timestep to write.)
		size_t count[] = {1,(size_t)(dims)};
		size_t start[] = {0,0};

		// Write the data to the file.
		if ( nc_put_vara_double(fid,id,start,count,data) != NC_NOERR ) return NETCDF_ERR;

		return NETCDF_OK;
	}
	int NetCDFPutVar(int fid, int id, int dims, float *data) {
		// These settings tell netcdf to write one timestep of data. (The
		// setting of start[0] inside the loop below tells netCDF which
		// timestep to write.)
		size_t count[] = {1,(size_t)(dims)};
		size_t start[] = {0,0};

		// Write the data to the file.
		if ( nc_put_vara_float(fid,id,start,count,data) != NC_NOERR ) return NETCDF_ERR;

		return NETCDF_OK;
	}
	int NetCDFPutVar(int fid, int id, int dims[], double *data) {
		// These settings tell netcdf to write one timestep of data. (The
		// setting of start[0] inside the loop below tells netCDF which
		// timestep to write.)
		size_t count[] = {1,(size_t)(dims[2]),(size_t)(dims[1]),(size_t)(dims[0])};
		size_t start[] = {0,0,0,0};

		// Write the data to the file.
		if ( nc_put_vara_double(fid,id,start,count,data) != NC_NOERR ) return NETCDF_ERR;

		return NETCDF_OK;
	}
	int NetCDFPutVar(int fid, int id, int dims[], float *data) {
		// These settings tell netcdf to write one timestep of data. (The
		// setting of start[0] inside the loop below tells netCDF which
		// timestep to write.)
		size_t count[] = {1,(size_t)(dims[2]),(size_t)(dims[1]),(size_t)(dims[0])};
		size_t start[] = {0,0,0,0};

		// Write the data to the file.
		if ( nc_put_vara_float(fid,id,start,count,data) != NC_NOERR ) return NETCDF_ERR;
			
		return NETCDF_OK;
	}
	int NetCDFGetVar(int fid, const char *varname, const int n, double *out, bool m2zero) {
		int varid, vartype;
		// Get the variable id based on its name
		if ( nc_inq_varid(fid,varname,&varid) != NC_NOERR ) return NETCDF_ERR;
		// Get the variable type
		if ( nc_inq_vartype(fid,varid,&vartype) != NC_NOERR ) return NETCDF_ERR;
		// Read the data according to its type
		switch (vartype) {
			case NC_FLOAT:
				{
					std::vector<float> aux(n);
					if ( nc_get_var_float(fid,varid,&aux[0]) != NC_NOERR ) return NETCDF_ERR;
					std::copy(&aux[0],&aux[0]+n,out);
				}
				break;
			case NC_DOUBLE:
				if ( nc_get_var_double(fid,varid,out) != NC_NOERR ) return NETCDF_ERR;
				break;
			default:
				return NETCDF_ERR;
				break;
		}
		// Deal with the missing value (convert it to zeros)
		if (m2zero) {
			#pragma omp parallel
			{
			for (int ii = OMP_THREAD_NUM; ii < n; ii += OMP_NUM_THREADS)
				out[ii] = (out[ii] >= MISSING_VALUE) ? 0. : out[ii];
			}			
		}
		return NETCDF_OK;
	}
	int NetCDFGetVar(int fid, const char *varname, const int n, float *out, bool m2zero) {
		int varid, vartype;
		// Get the variable id based on its name
		if ( nc_inq_varid(fid,varname,&varid) != NC_NOERR ) return NETCDF_ERR;
		// Get the variable type
		if ( nc_inq_vartype(fid,varid,&vartype) != NC_NOERR ) return NETCDF_ERR;
		// Read the data according to its type
		switch (vartype) {
			case NC_FLOAT:
				if ( nc_get_var_float(fid,varid,out) != NC_NOERR ) return NETCDF_ERR;
				break;
			case NC_DOUBLE:
				{
					std::vector<double> aux(n);
					if ( nc_get_var_double(fid,varid,&aux[0]) != NC_NOERR ) return NETCDF_ERR;
					std::copy(&aux[0],&aux[0]+n,out);
				}
				break;
			default:
				return NETCDF_ERR;
				break;
		}
		// Deal with the missing value (convert it to zeros)
		if (m2zero) {
			#pragma omp parallel
			{
			for (int ii = OMP_THREAD_NUM; ii < n; ii += OMP_NUM_THREADS)
				out[ii] = (out[ii] >= MISSING_VALUE) ? 0. : out[ii];
			}			
		}
		return NETCDF_OK;
	}

	/* READNETCDF

		Reads a NetCDF4 file given the name of the file, the name of the variable
		to be read and its size.
	*/
	int readNetCDF(const char *fname, const char *varname, const int n, double *out, bool m2zero) {

		int fid;

		// Open file for reading
		if ( nc_open(fname,NC_NOWRITE,&fid) != NC_NOERR ) return NETCDF_ERR;

		// Read the data
		if ( NetCDFGetVar(fid,varname,n,out,m2zero) != NETCDF_OK ) return NETCDF_ERR;

		// Close the file
		if ( nc_close(fid) != NC_NOERR) return NETCDF_ERR;

		return NETCDF_OK;
	}
	int readNetCDF(const char *fname, const char *varname, const int n, float *out, bool m2zero) {

		int fid;

		// Open file for reading
		if ( nc_open(fname,NC_NOWRITE,&fid) != NC_NOERR ) return NETCDF_ERR;

		// Read the data
		if ( NetCDFGetVar(fid,varname,n,out,m2zero) != NETCDF_OK ) return NETCDF_ERR;

		// Close the file
		if ( nc_close(fid) != NC_NOERR) return NETCDF_ERR;

		return NETCDF_OK;
	}
	int readNetCDF(const char *fname, const char *varname, std::vector<double> &out) {
		return readNetCDF(fname,varname,(int)(out.size()),&out[0],true);
	}
	int readNetCDF(const char *fname, const char *varname, std::vector<float> &out) {
		return readNetCDF(fname,varname,(int)(out.size()),&out[0],true);
	}
	int readNetCDF(const char *fname, const char *varname, field::Field<double> &f) {
		return readNetCDF(fname,varname,f.get_n(),f.data(),true);
	}
	int readNetCDF(const char *fname, const char *varname, field::Field<float> &f) {
		return readNetCDF(fname,varname,f.get_n(),f.data(),true);
	}
	int readNetCDF(const char *fname, std::string *varname, field::Field<double> &f) {

		int retval;

		// Define auxiliary variable
		std::vector<std::vector<double>> aux_var(f.get_m());

		// Load each variable in the auxiliary vector, do not convert empties to zeros
		for (int vid = 0; vid < f.get_m(); ++vid) {
			aux_var[vid].resize(f.get_n());
			retval = readNetCDF(fname,varname[vid].c_str(),f.get_n(),&aux_var[vid][0],false);
		}

		// Set output field and eliminate the missing variables
		#pragma omp parallel
		{
		for (int ii = OMP_THREAD_NUM; ii < f.get_n(); ii += OMP_NUM_THREADS)
			#pragma loop_count min(1), max(16), avg(3) // for vectorization
			for (int jj = 0; jj < f.get_m(); ++jj)
				f[ii][jj] = (aux_var[jj][ii] >= MISSING_VALUE) ? 0. : aux_var[jj][ii];
		}

		return retval;
	}
	int readNetCDF(const char *fname, std::string *varname, field::Field<float> &f) {

		int retval;

		// Define auxiliary variable
		std::vector<std::vector<float>> aux_var(f.get_m());

		// Load each variable in the auxiliary vector, do not convert empties to zeros
		for (int vid = 0; vid < f.get_m(); ++vid) {
			aux_var[vid].resize(f.get_n());
			retval = readNetCDF(fname,varname[vid].c_str(),f.get_n(),&aux_var[vid][0],false);
		}

		// Set output field and eliminate the missing variables
		#pragma omp parallel
		{
		for (int ii = OMP_THREAD_NUM; ii < f.get_n(); ii += OMP_NUM_THREADS)
			#pragma loop_count min(1), max(16), avg(3) // for vectorization
			for (int jj = 0; jj < f.get_m(); ++jj)
				f[ii][jj] = (aux_var[jj][ii] >= MISSING_VALUE) ? 0. : aux_var[jj][ii];
		}

		return retval;
	}

	/* WRITENETCDF

		Writes a NetCDF4 file given the variable, the name of the file, the name of the variable
		to be read and its size.
	*/
	int writeNetCDF(const char *fname, const char *varname, int dims[], double *lon, double *lat, 
		double *depth, double *data) {
		
		int fid, lon_id, lat_id, dep_id, tim_id;
		int varid, vlon_id, vlat_id, vdep_id;

		// Create the file
		if ( NetCDFCreate(fname,fid) != NETCDF_OK) return NETCDF_ERR;

		// Define the dimensions. NetCDF will hand back 2 IDs for each.
		if ( NetCDFDefXY(fid,lon_id,vlon_id,"x","nav_lon",NC_DOUBLE,dims[0],"Longitude","degrees_east")  != NETCDF_OK ) return NETCDF_ERR;
		if ( NetCDFDefXY(fid,lat_id,vlat_id,"y","nav_lat",NC_DOUBLE,dims[1],"Longitude","degrees_north") != NETCDF_OK ) return NETCDF_ERR;
		if ( NetCDFDefZ(fid,dep_id,vdep_id,"deptht","deptht",NC_DOUBLE,dims[2],"meter","down")           != NETCDF_OK ) return NETCDF_ERR;
		if ( NetCDFDefT(fid,tim_id,"time_counter")                                                       != NETCDF_OK ) return NETCDF_ERR;

		// Define the variables
		int dimvar[4] = {tim_id,dep_id,lat_id,lon_id};
		if ( NetCDFDefVar(fid,varid,varname,NC_DOUBLE,4,dimvar) != NETCDF_OK ) return NETCDF_ERR;
		
		// End define mode
		if ( nc_enddef(fid) != NC_NOERR ) return NETCDF_ERR;

		// Write the data to the file.
		if ( NetCDFPutDim(fid,vlon_id,&lon[0])     != NETCDF_OK ) return NETCDF_ERR;
		if ( NetCDFPutDim(fid,vlat_id,&lat[0])     != NETCDF_OK ) return NETCDF_ERR;
		if ( NetCDFPutDim(fid,vdep_id,&depth[0])   != NETCDF_OK ) return NETCDF_ERR;
		if ( NetCDFPutVar(fid,varid,dims,&data[0]) != NETCDF_OK ) return NETCDF_ERR;

		// Close the file
		if ( nc_close(fid) != NC_NOERR) return NETCDF_ERR;

		return NETCDF_OK;
	}
	int writeNetCDF(const char *fname, const char *varname, int dims[], float *lon, float *lat, 
		float *depth, float *data) {

		int fid, lon_id, lat_id, dep_id, tim_id;
		int varid, vlon_id, vlat_id, vdep_id;

		// Create the file
		if ( NetCDFCreate(fname,fid) != NETCDF_OK) return NETCDF_ERR;

		// Define the dimensions. NetCDF will hand back 2 IDs for each.
		if ( NetCDFDefXY(fid,lon_id,vlon_id,"x","nav_lon",NC_FLOAT,dims[0],"Longitude","degrees_east")  != NETCDF_OK ) return NETCDF_ERR;
		if ( NetCDFDefXY(fid,lat_id,vlat_id,"y","nav_lat",NC_FLOAT,dims[1],"Longitude","degrees_north") != NETCDF_OK ) return NETCDF_ERR;
		if ( NetCDFDefZ(fid,dep_id,vdep_id,"deptht","deptht",NC_FLOAT,dims[2],"meter","down")           != NETCDF_OK ) return NETCDF_ERR;
		if ( NetCDFDefT(fid,tim_id,"time_counter")                                                      != NETCDF_OK ) return NETCDF_ERR;

		// Define the variables
		int dimvar[4] = {tim_id,dep_id,lat_id,lon_id};
		if ( NetCDFDefVar(fid,varid,varname,NC_FLOAT,4,dimvar) != NETCDF_OK ) return NETCDF_ERR;
		
		// End define mode
		if ( nc_enddef(fid) != NC_NOERR ) return NETCDF_ERR;

		// Write the data to the file.
		if ( NetCDFPutDim(fid,vlon_id,&lon[0])     != NETCDF_OK ) return NETCDF_ERR;
		if ( NetCDFPutDim(fid,vlat_id,&lat[0])     != NETCDF_OK ) return NETCDF_ERR;
		if ( NetCDFPutDim(fid,vdep_id,&depth[0])   != NETCDF_OK ) return NETCDF_ERR;
		if ( NetCDFPutVar(fid,varid,dims,&data[0]) != NETCDF_OK ) return NETCDF_ERR;

		// Close the file
		if ( nc_close(fid) != NC_NOERR) return NETCDF_ERR;

		return NETCDF_OK;		
	}
	int writeNetCDF(const char *fname, std::string *varname, int nvars, int dims[], double *lon, double *lat, 
		double *depth, double **data) {

		int fid, lon_id, lat_id, dep_id, tim_id;
		int vlon_id, vlat_id, vdep_id;

		// Create the file
		if ( NetCDFCreate(fname,fid) != NETCDF_OK) return NETCDF_ERR;

		// Define the dimensions. NetCDF will hand back 2 IDs for each.
		if ( NetCDFDefXY(fid,lon_id,vlon_id,"x","nav_lon",NC_DOUBLE,dims[0],"Longitude","degrees_east")  != NETCDF_OK ) return NETCDF_ERR;
		if ( NetCDFDefXY(fid,lat_id,vlat_id,"y","nav_lat",NC_DOUBLE,dims[1],"Longitude","degrees_north") != NETCDF_OK ) return NETCDF_ERR;
		if ( NetCDFDefZ(fid,dep_id,vdep_id,"deptht","deptht",NC_DOUBLE,dims[2],"meter","down")           != NETCDF_OK ) return NETCDF_ERR;
		if ( NetCDFDefT(fid,tim_id,"time_counter")                                                       != NETCDF_OK ) return NETCDF_ERR;

		// Define the variables
		int dimvar[4] = {tim_id,dep_id,lat_id,lon_id};
		std::vector<int> varid(nvars);
		for (int vid = 0; vid < nvars; ++vid)
			if ( NetCDFDefVar(fid,varid[vid],varname[vid].c_str(),NC_DOUBLE,4,dimvar) != NETCDF_OK ) return NETCDF_ERR;
		
		// End define mode
		if ( nc_enddef(fid) != NC_NOERR ) return NETCDF_ERR;

		// Write the data to the file.
		if ( NetCDFPutDim(fid,vlon_id,&lon[0])     != NETCDF_OK ) return NETCDF_ERR;
		if ( NetCDFPutDim(fid,vlat_id,&lat[0])     != NETCDF_OK ) return NETCDF_ERR;
		if ( NetCDFPutDim(fid,vdep_id,&depth[0])   != NETCDF_OK ) return NETCDF_ERR;
		for (int vid = 0; vid < nvars; ++vid)
			if ( NetCDFPutVar(fid,varid[vid],dims,&data[vid][0]) != NETCDF_OK ) return NETCDF_ERR;

		// Close the file
		if ( nc_close(fid) != NC_NOERR) return NETCDF_ERR;

		return NETCDF_OK;		
	}
	int writeNetCDF(const char *fname, std::string *varname, int nvars, int dims[], float *lon, float *lat, 
		float *depth, float **data) {

		int fid, lon_id, lat_id, dep_id, tim_id;
		int vlon_id, vlat_id, vdep_id;

		// Create the file
		if ( NetCDFCreate(fname,fid) != NETCDF_OK) return NETCDF_ERR;

		// Define the dimensions. NetCDF will hand back 2 IDs for each.
		if ( NetCDFDefXY(fid,lon_id,vlon_id,"x","nav_lon",NC_FLOAT,dims[0],"Longitude","degrees_east")  != NETCDF_OK ) return NETCDF_ERR;
		if ( NetCDFDefXY(fid,lat_id,vlat_id,"y","nav_lat",NC_FLOAT,dims[1],"Longitude","degrees_north") != NETCDF_OK ) return NETCDF_ERR;
		if ( NetCDFDefZ(fid,dep_id,vdep_id,"deptht","deptht",NC_FLOAT,dims[2],"meter","down")           != NETCDF_OK ) return NETCDF_ERR;
		if ( NetCDFDefT(fid,tim_id,"time_counter")                                                      != NETCDF_OK ) return NETCDF_ERR;

		// Define the variables
		int dimvar[4] = {tim_id,dep_id,lat_id,lon_id};
		std::vector<int> varid(nvars);
		for (int vid = 0; vid < nvars; ++vid)
			if ( NetCDFDefVar(fid,varid[vid],varname[vid].c_str(),NC_FLOAT,4,dimvar) != NETCDF_OK ) return NETCDF_ERR;
		
		// End define mode
		if ( nc_enddef(fid) != NC_NOERR ) return NETCDF_ERR;

		// Write the data to the file.
		if ( NetCDFPutDim(fid,vlon_id,&lon[0])     != NETCDF_OK ) return NETCDF_ERR;
		if ( NetCDFPutDim(fid,vlat_id,&lat[0])     != NETCDF_OK ) return NETCDF_ERR;
		if ( NetCDFPutDim(fid,vdep_id,&depth[0])   != NETCDF_OK ) return NETCDF_ERR;
		for (int vid = 0; vid < nvars; ++vid)
			if ( NetCDFPutVar(fid,varid[vid],dims,&data[vid][0]) != NETCDF_OK ) return NETCDF_ERR;

		// Close the file
		if ( nc_close(fid) != NC_NOERR) return NETCDF_ERR;

		return NETCDF_OK;		
	}

	int writeNetCDF(const char *fname, const char *varname, int dims[], double *lon, double *lat, 
		double *depth, field::Field<double> &f) {
		// Write NetCDF
		int retval;
		if (f.get_m() == 1) {
			// We can use the writeNetCDF API to write a single file containing a variable
			retval = writeNetCDF(fname,varname,dims,&lon[0],&lat[0],&depth[0],&f[0][0]);
		} else {
			// For variables with multiple components, we need to create variable names
			std::vector<std::string> newvarname(f.get_m());
			std::vector<double*> ff(f.get_m());
			// Generate variable name and array
			for (int mm = 0; mm < f.get_m(); ++mm) {
				newvarname[mm] = std::string(varname) + std::string("_") + std::to_string(mm);
				std::replace(newvarname[mm].begin(),newvarname[mm].end(),' ','_'); // Replace white spaces
				// Fill the array
				ff[mm] = new double[f.get_n()];
				for (int nn = 0; nn < f.get_n(); ++nn)
					ff[mm][nn] = f[nn][mm];
			}
			// Write NetCDF
			retval = writeNetCDF(fname,newvarname.data(),f.get_m(),dims,&lon[0],&lat[0],&depth[0],ff.data());
			// Deallocate
			for (int mm = 0; mm < f.get_m(); ++mm)
				delete [] ff[mm];
		}

		return retval;
	}
	int writeNetCDF(const char *fname, const char *varname, int dims[], float *lon, float *lat, 
		float *depth, field::Field<float> &f) {
		// Write NetCDF
		int retval;
		if (f.get_m() == 1) {
			// We can use the writeNetCDF API to write a single file containing a variable
			retval = writeNetCDF(fname,varname,dims,&lon[0],&lat[0],&depth[0],&f[0][0]);
		} else {
			// For variables with multiple components, we need to create variable names
			std::vector<std::string> newvarname(f.get_m());
			std::vector<float*> ff(f.get_m());
			// Generate variable name and array
			for (int mm = 0; mm < f.get_m(); ++mm) {
				newvarname[mm] = std::string(varname) + std::string("_") + std::to_string(mm);
				std::replace(newvarname[mm].begin(),newvarname[mm].end(),' ','_'); // Replace white spaces
				// Fill the array
				ff[mm] = new float[f.get_n()];
				for (int nn = 0; nn < f.get_n(); ++nn)
					ff[mm][nn] = f[nn][mm];
			}
			// Write NetCDF
			retval = writeNetCDF(fname,newvarname.data(),f.get_m(),dims,&lon[0],&lat[0],&depth[0],ff.data());
			// Deallocate
			for (int mm = 0; mm < f.get_m(); ++mm)
				delete [] ff[mm];
		}

		return retval;
	}

	/* WRITENETCDFPROFILE

		Writes a NetCDF4 file given the variable, the name of the file, the name of the variable
		to be read and its size.
	*/
	int writeNetCDFProfile(const char *fname, const char *varname, int dims, double *lon, double *lat, 
		double *depth, double *data) {
		
		int fid, loid, laid, deid, tid;
		int varid, vloid, vlaid, vdeid;

		// Create the file
		if ( NetCDFCreate(fname,fid) != NETCDF_OK) return NETCDF_ERR;

		// Define the dimensions. NetCDF will hand back 2 IDs for each.
		if ( NetCDFDefZ(fid,deid,vdeid,"deptht","deptht",NC_DOUBLE,dims,"meter","down")           != NETCDF_OK ) return NETCDF_ERR;
		if ( NetCDFDefXY(fid,loid,vloid,"x","nav_lon",NC_DOUBLE,dims,"Longitude","degrees_east")  != NETCDF_OK ) return NETCDF_ERR;
		if ( NetCDFDefXY(fid,laid,vlaid,"y","nav_lat",NC_DOUBLE,dims,"Longitude","degrees_north") != NETCDF_OK ) return NETCDF_ERR;		
		if ( NetCDFDefT(fid,tid,"time_counter")                                                   != NETCDF_OK ) return NETCDF_ERR;

		// Define the variables
		int dimvar[2] = {tid,vdeid};
		if ( NetCDFDefVar(fid,varid,varname,NC_DOUBLE,2,dimvar)   != NETCDF_OK ) return NETCDF_ERR;
		
		// End define mode
		if ( nc_enddef(fid) != NC_NOERR ) return NETCDF_ERR;

		// Write the data to the file.
		if ( NetCDFPutDim(fid,vloid,&lon[0])       != NETCDF_OK ) return NETCDF_ERR;
		if ( NetCDFPutDim(fid,vlaid,&lat[0])       != NETCDF_OK ) return NETCDF_ERR;
		if ( NetCDFPutDim(fid,vdeid,&depth[0])     != NETCDF_OK ) return NETCDF_ERR;
		if ( NetCDFPutVar(fid,varid,dims,&data[0]) != NETCDF_OK ) return NETCDF_ERR;

		// Close the file
		if ( nc_close(fid) != NC_NOERR) return NETCDF_ERR;

		return NETCDF_OK;
	}
	int writeNetCDFProfile(const char *fname, const char *varname, int dims, float *lon, float *lat, 
		float *depth, float *data) {
		
		int fid, loid, laid, deid, tid;
		int varid, vloid, vlaid, vdeid;

		// Create the file
		if ( NetCDFCreate(fname,fid) != NETCDF_OK) return NETCDF_ERR;

		// Define the dimensions. NetCDF will hand back 2 IDs for each.
		if ( NetCDFDefZ(fid,deid,vdeid,"deptht","deptht",NC_FLOAT,dims,"meter","down") != NETCDF_OK ) return NETCDF_ERR;
		if ( NetCDFDefXY(fid,loid,vloid,"x","nav_lon",NC_FLOAT,dims,"Longitude","degrees_east")  != NETCDF_OK ) return NETCDF_ERR;
		if ( NetCDFDefXY(fid,laid,vlaid,"y","nav_lat",NC_FLOAT,dims,"Longitude","degrees_north") != NETCDF_OK ) return NETCDF_ERR;				
		if ( NetCDFDefT(fid,tid,"time_counter")                                       != NETCDF_OK ) return NETCDF_ERR;

		// Define the variables
		int dimvar[2] = {tid,vdeid};
		if ( NetCDFDefVar(fid,varid,varname,NC_FLOAT,2,dimvar)   != NETCDF_OK ) return NETCDF_ERR;
		
		// End define mode
		if ( nc_enddef(fid) != NC_NOERR ) return NETCDF_ERR;

		// Write the data to the file.
		if ( NetCDFPutDim(fid,vloid,&lon[0])       != NETCDF_OK ) return NETCDF_ERR;
		if ( NetCDFPutDim(fid,vlaid,&lat[0])       != NETCDF_OK ) return NETCDF_ERR;
		if ( NetCDFPutDim(fid,vdeid,&depth[0])     != NETCDF_OK ) return NETCDF_ERR;
		if ( NetCDFPutVar(fid,varid,dims,&data[0]) != NETCDF_OK ) return NETCDF_ERR;

		// Close the file
		if ( nc_close(fid) != NC_NOERR) return NETCDF_ERR;

		return NETCDF_OK;
	}
	int writeNetCDFProfile(const char *fname, std::string *varname, int nvars, int dims, double *lon, double *lat, 
		double *depth, double **data) {
		
		int fid, loid, laid, deid, tid;
		int vloid, vlaid, vdeid;

		// Create the file
		if ( NetCDFCreate(fname,fid) != NETCDF_OK) return NETCDF_ERR;

		// Define the dimensions. NetCDF will hand back 2 IDs for each.
		if ( NetCDFDefZ(fid,deid,vdeid,"deptht","deptht",NC_DOUBLE,dims,"meter","down")          != NETCDF_OK ) return NETCDF_ERR;
		if ( NetCDFDefXY(fid,loid,vloid,"x","nav_lon",NC_FLOAT,dims,"Longitude","degrees_east")  != NETCDF_OK ) return NETCDF_ERR;
		if ( NetCDFDefXY(fid,laid,vlaid,"y","nav_lat",NC_FLOAT,dims,"Longitude","degrees_north") != NETCDF_OK ) return NETCDF_ERR;				
		if ( NetCDFDefT(fid,tid,"time_counter")                                                  != NETCDF_OK ) return NETCDF_ERR;

		// Define the variables
		int dimvar[2] = {tid,vdeid};
		std::vector<int> varid(nvars);
		for (int vid = 0; vid < nvars; ++vid)
			if ( NetCDFDefVar(fid,varid[vid],varname[vid].c_str(),NC_DOUBLE,2,dimvar) != NETCDF_OK ) return NETCDF_ERR;
		
		// End define mode
		if ( nc_enddef(fid) != NC_NOERR ) return NETCDF_ERR;

		// Write the data to the file.
		if ( NetCDFPutDim(fid,vloid,&lon[0])   != NETCDF_OK ) return NETCDF_ERR;
		if ( NetCDFPutDim(fid,vlaid,&lat[0])   != NETCDF_OK ) return NETCDF_ERR;
		if ( NetCDFPutDim(fid,vdeid,&depth[0]) != NETCDF_OK ) return NETCDF_ERR;
		for (int vid = 0; vid < nvars; ++vid)
			if ( NetCDFPutVar(fid,varid[vid],dims,&data[vid][0]) != NETCDF_OK ) return NETCDF_ERR;

		// Close the file
		if ( nc_close(fid) != NC_NOERR) return NETCDF_ERR;

		return NETCDF_OK;
	}
	int writeNetCDFProfile(const char *fname, std::string *varname, int nvars, int dims, float *lon, float *lat, 
		float *depth, float **data) {
		
		int fid, loid, laid, deid, tid;
		int vloid, vlaid, vdeid;

		// Create the file
		if ( NetCDFCreate(fname,fid) != NETCDF_OK) return NETCDF_ERR;

		// Define the dimensions. NetCDF will hand back 2 IDs for each.
		if ( NetCDFDefZ(fid,deid,vdeid,"deptht","deptht",NC_FLOAT,dims,"meter","down") != NETCDF_OK ) return NETCDF_ERR;
		if ( NetCDFDefXY(fid,loid,vloid,"x","nav_lon",NC_FLOAT,dims,"Longitude","degrees_east")  != NETCDF_OK ) return NETCDF_ERR;
		if ( NetCDFDefXY(fid,laid,vlaid,"y","nav_lat",NC_FLOAT,dims,"Longitude","degrees_north") != NETCDF_OK ) return NETCDF_ERR;				
		if ( NetCDFDefT(fid,tid,"time_counter")                                        != NETCDF_OK ) return NETCDF_ERR;

		// Define the variables
		int dimvar[2] = {tid,vdeid};
		std::vector<int> varid(nvars);
		for (int vid = 0; vid < nvars; ++vid)
			if ( NetCDFDefVar(fid,varid[vid],varname[vid].c_str(),NC_FLOAT,2,dimvar) != NETCDF_OK ) return NETCDF_ERR;
		
		// End define mode
		if ( nc_enddef(fid) != NC_NOERR ) return NETCDF_ERR;

		// Write the data to the file.
		if ( NetCDFPutDim(fid,vloid,&lon[0])   != NETCDF_OK ) return NETCDF_ERR;
		if ( NetCDFPutDim(fid,vlaid,&lat[0])   != NETCDF_OK ) return NETCDF_ERR;
		if ( NetCDFPutDim(fid,vdeid,&depth[0]) != NETCDF_OK ) return NETCDF_ERR;
		for (int vid = 0; vid < nvars; ++vid)
			if ( NetCDFPutVar(fid,varid[vid],dims,&data[vid][0]) != NETCDF_OK ) return NETCDF_ERR;

		// Close the file
		if ( nc_close(fid) != NC_NOERR) return NETCDF_ERR;

		return NETCDF_OK;
	}
	int writeNetCDFProfile(const char *fname, const char *varname, int dims, double *lon, double *lat, 
		double *depth, field::Field<double> &f) {
		// Write NetCDF
		int retval;
		if (f.get_m() == 1) {
			// We can use the writeNetCDF API to write a single file containing a variable
			retval = writeNetCDFProfile(fname,varname,dims,&lon[0],&lat[0],&depth[0],&f[0][0]);
		} else {
			// For variables with multiple components, we need to create variable names
			std::vector<std::string> newvarname(f.get_m());
			std::vector<double*> ff(f.get_m());
			// Generate variable name and array
			for (int mm = 0; mm < f.get_m(); ++mm) {
				newvarname[mm] = std::string(varname) + std::string("_") + std::to_string(mm);
				std::replace(newvarname[mm].begin(),newvarname[mm].end(),' ','_'); // Replace white spaces
				// Fill the array
				ff[mm] = new double[f.get_n()];
				for (int nn = 0; nn < f.get_n(); ++nn)
					ff[mm][nn] = f[nn][mm];
			}
			// Write NetCDF
			retval = writeNetCDFProfile(fname,newvarname.data(),f.get_m(),dims,&lon[0],&lat[0],&depth[0],ff.data());
			// Deallocate
			for (int mm = 0; mm < f.get_m(); ++mm)
				delete [] ff[mm];
		}
		return retval;
	}
	int writeNetCDFProfile(const char *fname, const char *varname, int dims, float *lon, float *lat, 
		float *depth, field::Field<float> &f) {
		// Write NetCDF
		int retval;
		if (f.get_m() == 1) {
			// We can use the writeNetCDF API to write a single file containing a variable
			retval = writeNetCDFProfile(fname,varname,dims,&lon[0],&lat[0],&depth[0],&f[0][0]);
		} else {
			// For variables with multiple components, we need to create variable names
			std::vector<std::string> newvarname(f.get_m());
			std::vector<float*> ff(f.get_m());
			// Generate variable name and array
			for (int mm = 0; mm < f.get_m(); ++mm) {
				newvarname[mm] = std::string(varname) + std::string("_") + std::to_string(mm);
				std::replace(newvarname[mm].begin(),newvarname[mm].end(),' ','_'); // Replace white spaces
				// Fill the array
				ff[mm] = new float[f.get_n()];
				for (int nn = 0; nn < f.get_n(); ++nn)
					ff[mm][nn] = f[nn][mm];
			}
			// Write NetCDF
			retval = writeNetCDFProfile(fname,newvarname.data(),f.get_m(),dims,&lon[0],&lat[0],&depth[0],ff.data());
			// Deallocate
			for (int mm = 0; mm < f.get_m(); ++mm)
				delete [] ff[mm];
		}
		return retval;
	}
}