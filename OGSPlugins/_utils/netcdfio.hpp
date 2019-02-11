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

#include "field.h"

namespace NetCDF
{
	/* Read NetCDF Routines */
	double              *readNetCDF(const char *fname, const char *varname, const int n);
	field::Field<double> readNetCDF2F(const char *fname, const char *varname, const int n);
  field::Field<double> readNetCDF2F3(const char *fname, const char *vname1, const char *vname2, const char *vname3, const int n);
}

#endif