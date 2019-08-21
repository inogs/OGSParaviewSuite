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

#ifndef NETCDFIO_H
#define NETCDFIO_H

#include "V3.h"
#include "field.h"

#include <vector>
#include <string>

#define NETCDF_ERR 0
#define NETCDF_OK  1
#define MISSING_VALUE 1.e20

namespace NetCDF
{
	/* Read NetCDF Routines */
  int readNetCDF(const char *fname, const char *varname, const int n, double *out, bool m2zero);
  int readNetCDF(const char *fname, const char *varname, const int n, float *out, bool m2zero);
  int readNetCDF(const char *fname, const char *varname, std::vector<double> &out);
  int readNetCDF(const char *fname, const char *varname, std::vector<float> &out);
  int readNetCDF(const char *fname, const char *varname, field::Field<double> &f);
  int readNetCDF(const char *fname, const char *varname, field::Field<float> &f);
  int readNetCDF(const char *fname, std::string *varname, field::Field<double> &f);
  int readNetCDF(const char *fname, std::string *varname, field::Field<float> &f);
  /* Write NetCDF Routines */
  int writeNetCDF(const char *fname, const char *varname, int dims[], double *lon, double *lat, double *depth, double *data);
  int writeNetCDF(const char *fname, const char *varname, int dims[], float *lon, float *lat, float *depth, float *data);
  int writeNetCDF(const char *fname, std::string *varname, int nvars, int dims[], double *lon, double *lat, double *depth, double **data);
  int writeNetCDF(const char *fname, std::string *varname, int nvars, int dims[], float *lon, float *lat, float *depth, float **data);
  int writeNetCDF(const char *fname, std::string *varname, int nvars, int dims[], v3::V3v &xyz, double **data);
  int writeNetCDF(const char *fname, std::string *varname, int nvars, int dims[], v3::V3v &xyz, float **data);
  int writeNetCDF(const char *fname, const char *varname, int dims[], v3::V3v &xyz, field::Field<double> &f);
  int writeNetCDF(const char *fname, const char *varname, int dims[], v3::V3v &xyz, field::Field<float> &f);
  int writeNetCDFProfile(const char *fname, const char *varname, int dims, double *lon, double *lat, double *depth, double *data);
  int writeNetCDFProfile(const char *fname, const char *varname, int dims, float *lon, float *lat, float *depth, float *data);
  int writeNetCDFProfile(const char *fname, std::string *varname, int nvars, int dims, double *lon, double *lat, double *depth, double **data);
  int writeNetCDFProfile(const char *fname, std::string *varname, int nvars, int dims, float *lon, float *lat, float *depth, float **data);
  int writeNetCDFProfile(const char *fname, std::string *varname, int nvars, int dims, v3::V3v &xyz, double **data);
  int writeNetCDFProfile(const char *fname, std::string *varname, int nvars, int dims, v3::V3v &xyz, float **data);
  int writeNetCDFProfile(const char *fname, const char *varname, int dims, v3::V3v &xyz, field::Field<double> &f);
  int writeNetCDFProfile(const char *fname, const char *varname, int dims, v3::V3v &xyz, field::Field<float> &f);
}

#endif