/*=========================================================================

  Module:    OGS Mesh functions

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef OGSmesh_h
#define OGSmesh_h

/* Function prototypes */

extern "C" void writeOGSMesh(const char *fname, int nLon, int nLat, int nLev,
	double *Lon2Meters, double *Lat2Meters, double *nav_lev,
	double* basins_mask, double *coast_mask);

extern "C" void readOGSMesh(const char *fname, int *nLon, int *nLat, int *nLev,
	double **Lon2Meters, double **Lat2Meters, double **nav_lev,
	double **basins_mask, double **coast_mask);

#endif