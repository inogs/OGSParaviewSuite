/*=========================================================================

  Program:   Utilities
  Module:    vtkOperations.h

  Perform operations on vtk arrays that might be needed accross plugins.

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef VTKOPERATIONS_H
#define VTKOPERATIONS_H

/*
	These macros return the index needed to iterate an array
*/
#define VTKIND(ii,jj,kk,nx,ny)          ( (nx-1)*(ny-1)*(kk) + (nx-1)*(jj) + (ii) )
#define ARRIND(ii,jj,kk,nx,ny)          ( (nx)*(ny)*(kk) + (nx)*(jj) + (ii) )
#define STTIND(bId,cId,kk,sId,ns,nz,nc) ( ns*nz*nc*bId + ns*nz*cId + ns*kk + sId )


/* CREATERECTILINEARGRID

	Creates a VTK rectilinear grid given the mesh dimensions and the 
	conversion to meters.
*/
void createRectilinearGrid(int nLon, int nLat, int nLev,
	double *Lon2Meters, double* Lat2Meters, double *nav_lev, double dps,
	vtkRectilinearGrid *rgrid);

/* CREATEVTKSCAF
	
	Creates a vtkFloatArray for a scalar field given the array name, its dimensions
	and an additional array to fill.
*/
vtkFloatArray *createVTKscaf(const char *name, int nx, int ny, int nz, double *array);
vtkFloatArray *createVTKscaf(const char *name, int n, double *array);
vtkFloatArray *createVTKscaffrom2d(const char *name, int nx, int ny, int nz, double *array);

/* CREATEVTKVECF3
	
	Creates a vtkFloatArray for a 3 component vector field given the array name, 
	its dimensions and an additional array to fill.
*/
vtkFloatArray *createVTKvecf3(const char *name, int nx, int ny, int nz, 
	double *a1, double *a2, double *a3);
vtkFloatArray *createVTKvecf3(const char *name, int n, double *a1, double *a2, double *a3);

/* CREATEVTKSTRF

	Creates a vtk string array given the name 
	and data to be converted.
*/
vtkStringArray *createVTKstrf(const char *name,const char *data);

/* GET COORDINATES

	Computes the cell/point coordinates given a VTKDataSet
*/
vtkFloatArray *getCellCoordinates(const char *varname, vtkDataSet *mesh);
vtkFloatArray *getPointCoordinates(const char *varname, vtkDataSet *mesh);

#endif