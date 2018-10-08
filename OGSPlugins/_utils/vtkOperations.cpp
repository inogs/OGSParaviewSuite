/*=========================================================================

  Program:   Utilities
  Module:    vtkOperations.cpp

  Perform operations on vtk arrays that might be needed accross plugins.

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkFloatArray.h"
#include "vtkStringArray.h"
#include "vtkRectilinearGrid.h"

/*
	This macro returns the index needed to iterate an array
*/
#define VTKIND(ii,jj,kk,nx,ny) ( (nx-1)*(ny-1)*(kk) + (nx-1)*(jj) + (ii) )
#define ARRIND(ii,jj,kk,nx,ny) ( (nx)*(ny)*(kk) + (nx)*(jj) + (ii) )

/* CREATERECTILINEARGRID

	Creates a VTK rectilinear grid given the mesh dimensions and the 
	conversion to meters.
*/
void createRectilinearGrid(int nLon, int nLat, int nLev,
	double *Lon2Meters, double* Lat2Meters, double *nav_lev, double dps,
	vtkRectilinearGrid *rgrid) {

	// Set vtkFloatArrays
	vtkFloatArray *daLon = vtkFloatArray::New();
	for (int ii = 0; ii < nLon; ii++) daLon->InsertNextValue(1.*Lon2Meters[ii]);

	vtkFloatArray *daLat = vtkFloatArray::New();
	for (int ii = 0; ii < nLat; ii++) daLat->InsertNextValue(1.*Lat2Meters[ii]);

	vtkFloatArray *daLev = vtkFloatArray::New();
	for (int ii=0; ii<nLev; ii++) daLev->InsertNextValue(-dps*nav_lev[ii]);

	// Set rectilinear grid
	rgrid->SetDimensions(nLon,nLat,nLev);
	rgrid->SetXCoordinates(daLon);
	rgrid->SetYCoordinates(daLat);
	rgrid->SetZCoordinates(daLev);
}

/* CREATEVTKSCAF
	
	Creates a vtkFloatArray for a scalar field given the array name, its dimensions
	in x, y, z coordinates and an additional array to fill.
*/
vtkFloatArray *createVTKscaf(const char *name, int nx, int ny, int nz, double *array) {
	/*
		Create a vtk scalar field given the name, dimensions and
		the array to be converted.
	*/
	// Create float array
	vtkFloatArray *vtkArray = vtkFloatArray::New();
	vtkArray->SetName(name);
	vtkArray->SetNumberOfComponents(1); // Scalar field
	vtkArray->SetNumberOfTuples(nx*ny*nz);
	vtkArray->Fill(0.); // Preallocate vtkArray to zero

	// Fill the vtkArray with the values of the array
	if (array != NULL) {
		for (int kk = 0; kk < nz-1; kk++) {
			for (int jj = 0; jj < ny-1; jj++) {
				for (int ii = 0; ii < nx-1; ii++) {
					vtkArray->SetTuple1(VTKIND(ii,jj,kk,nx,ny),
						array[ARRIND(ii,jj,kk,nx,ny)]);
				}
			}
		}
	}

	return vtkArray;
}

/* CREATEVTKSCAF
	
	Creates a vtkFloatArray for a scalar field given the array name, its dimensions
	and an additional array to fill.
*/
vtkFloatArray *createVTKscaf(const char *name, int n, double *array) {
	/*
		Create a vtk scalar field given the name, dimensions and
		the array to be converted.
	*/
	// Create float array
	vtkFloatArray *vtkArray = vtkFloatArray::New();
	vtkArray->SetName(name);
	vtkArray->SetNumberOfComponents(1); // Scalar field
	vtkArray->SetNumberOfTuples(n);
	vtkArray->Fill(0.); // Preallocate vtkArray to zero

	// Fill the vtkArray with the values of the array
	if (array != NULL) {
		for (int ii = 0; ii < n; ii++)
			vtkArray->SetTuple1(ii,array[ii]);
	}

	return vtkArray;
}

/* CREATEVTKSCAF
	
	Creates a vtkFloatArray for a 3 component vector field given the array name, 
	its dimensions in x, y, z coordinates and an additional array to fill.
*/
vtkFloatArray *createVTKvecf3(const char *name, int nx, int ny, int nz, double *a1, double *a2, double *a3) {

	// Create float array
	vtkFloatArray *vtkArray = vtkFloatArray::New();
	vtkArray->SetName(name);
	vtkArray->SetNumberOfComponents(3); // Vectorial field
	vtkArray->SetNumberOfTuples(nx*ny*nz);
	vtkArray->Fill(0.); // Preallocate vtkArray to zero

	// Fill the vtkArray with the values of the array
	if (a1 != NULL) {
		for (int kk = 0; kk < nz-1; kk++) {
			for (int jj = 0; jj < ny-1; jj++) {
				for (int ii = 0; ii < nx-1; ii++) {
					vtkArray->SetTuple3(VTKIND(ii,jj,kk,nx,ny),
						a1[ARRIND(ii,jj,kk,nx,ny)],
						a2[ARRIND(ii,jj,kk,nx,ny)],
						a3[ARRIND(ii,jj,kk,nx,ny)]);
				}
			}
		}
	}

	return vtkArray;
}

/* CREATEVTKSCAF
	
	Creates a vtkFloatArray for a 3 component vector field given the array name, 
	its dimensions and an additional array to fill.
*/
vtkFloatArray *createVTKvecf3(const char *name, int n, double *a1, double *a2, double *a3) {

	// Create float array
	vtkFloatArray *vtkArray = vtkFloatArray::New();
	vtkArray->SetName(name);
	vtkArray->SetNumberOfComponents(3); // Vectorial field
	vtkArray->SetNumberOfTuples(n);
	vtkArray->Fill(0.); // Preallocate vtkArray to zero

	// Fill the vtkArray with the values of the array
	if (a1 != NULL) {
		for (int ii = 0; ii < n; ii++)
			vtkArray->SetTuple3(ii,a1[ii],a2[ii],a3[ii]);
	}

	return vtkArray;
}

/* CREATEVTKSTRF

	Creates a vtk string array given the name 
	and data to be converted.
*/
vtkStringArray *createVTKstrf(const char *name,const char *data) {
	// Create string array
	vtkStringArray *vtkArray = vtkStringArray::New();
	vtkArray->SetName(name);
	vtkArray->SetNumberOfTuples(1);

	// Set the value
	vtkArray->SetValue(0,data);

	return vtkArray;
}