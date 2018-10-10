/*=========================================================================

  Module:    OGS Mesh functions

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include <cstdlib>
#include <cstdio>
#include <string>
#include <cstring>

#include "OGSdefs.h"
#include "OGSmesh.h"

using namespace std;

/* OGS Mesh writer 

	Generates a binary file (.ogsmsh) containing the transformation of 
	longitude and latitude to meters and the dept. The file is structured
	as follows:
		nLon, nLat, nLev (int)    --> dimensions
		Lon2Meters       (double) --> conversion of longitude to meters
		Lat2Meters       (double) --> conversion of latitude to meters
		nav_lev          (double) --> depth
		basins_mask      (double) --> sub basins mask
		coast_mask		 (double) --> coast mask
*/
extern "C" void writeOGSMesh(const char *fname, int nLon, int nLat, int nLev,
	double *Lon2Meters, double *Lat2Meters, double *nav_lev,
	double *basins_mask, double *coast_mask) {
	// Open file for writing
	FILE *myfile;
	myfile = fopen(fname,"wb"); if (myfile == NULL) ERROR("Cannot open file.")
	// Write the dimensions
	if (fwrite(&nLon,INTSZ,1,myfile) != 1) ERROR("Error writing file") 
	if (fwrite(&nLat,INTSZ,1,myfile) != 1) ERROR("Error writing file") 
	if (fwrite(&nLev,INTSZ,1,myfile) != 1) ERROR("Error writing file")
	// Write the vectors
	if (fwrite(Lon2Meters,DBLSZ,nLon,myfile) != nLon) ERROR("Error writing file")
	if (fwrite(Lat2Meters,DBLSZ,nLat,myfile) != nLat) ERROR("Error writing file")
	if (fwrite(nav_lev,DBLSZ,nLev,myfile)    != nLev) ERROR("Error writing file")
	// Write the masks
	int ncells = nLon*nLat*nLev;
	if (fwrite(basins_mask,DBLSZ,ncells,myfile) != ncells) ERROR("Error writing file")
	if (fwrite(coast_mask,DBLSZ,ncells,myfile)  != ncells) ERROR("Error writing file")
	// Close the file
	fclose(myfile);
}

/* OGS Mesh reader

	Reads a binary file (.ogsmsh) containing the transformation of 
	longitude and latitude to meters and the dept. The file is structured
	as follows:
		nLon, nLat, nLev (int)    --> dimensions
		Lon2Meters       (double) --> conversion of longitude to meters
		Lat2Meters       (double) --> conversion of latitude to meters
		nav_lev          (double) --> depth 
		basins_mask      (double) --> sub basins mask
		coast_mask		 (double) --> coast mask
*/
extern "C" void readOGSMesh(const char *fname, int *nLon, int *nLat, int *nLev,
	double **Lon2Meters, double **Lat2Meters, double **nav_lev,
	double **basins_mask, double **coast_mask) {
	// Open file for writing
	FILE *myfile;
	myfile = fopen(fname,"rb"); if (myfile == NULL) ERROR("Cannot open file.")
	// Read the dimensions
	if (fread(nLon,INTSZ,1,myfile) != 1) ERROR("Error reading file") 
	if (fread(nLat,INTSZ,1,myfile) != 1) ERROR("Error reading file") 
	if (fread(nLev,INTSZ,1,myfile) != 1) ERROR("Error reading file")
	// Read the vectors
	*Lon2Meters = (double*)malloc((*nLon)*DBLSZ); if (*Lon2Meters == NULL) ERROR("Cannot allocate memory")
	*Lat2Meters = (double*)malloc((*nLat)*DBLSZ); if (*Lat2Meters == NULL) ERROR("Cannot allocate memory")
	*nav_lev    = (double*)malloc((*nLev)*DBLSZ); if (*nav_lev    == NULL) ERROR("Cannot allocate memory")
	if (fread(*Lon2Meters,DBLSZ,*nLon,myfile) != (*nLon)) ERROR("Error reading file")
	if (fread(*Lat2Meters,DBLSZ,*nLat,myfile) != (*nLat)) ERROR("Error reading file")
	if (fread(*nav_lev,DBLSZ,*nLev,myfile)    != (*nLev)) ERROR("Error reading file")
	// Read the masks
	int ncells = (*nLon)*(*nLat)*(*nLev); 
	*basins_mask = (double*)malloc(ncells*DBLSZ); if (*basins_mask == NULL) ERROR("Cannot allocate memory")
	*coast_mask  = (double*)malloc(ncells*DBLSZ); if (*coast_mask  == NULL) ERROR("Cannot allocate memory")
	if (fread(*basins_mask,DBLSZ,ncells,myfile) != ncells) ERROR("Error reading file")
	if (fread(*coast_mask,DBLSZ,ncells,myfile)  != ncells) ERROR("Error reading file")
	// Close the file
	fclose(myfile);
}