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

#include "vtkCell.h"
#include "vtkPoints.h"
#include "vtkCellData.h"
#include "vtkPointData.h"
#include "vtkFieldData.h"

#include "vtkFloatArray.h"
#include "vtkStringArray.h"

#include "vtkDataSet.h"
#include "vtkRectilinearGrid.h"

#include <vector>

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

/* CREATEVTKSCAFFROM2D
	
	Creates a vtkFloatArray for a scalar field given the array name, its dimensions
	in x, y, z coordinates and an additional 2D array in z to fill.
*/
vtkFloatArray *createVTKscaffrom2d(const char *name, int nx, int ny, int nz, double *array) {
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
						array[ARRIND(ii,jj,0,nx,ny)]);
				}
			}
		}
	}

	return vtkArray;
}

/* CREATEVTKVECF3
	
	Creates a vtkFloatArray for a 3 component vector field given the array name, 
	its dimensions in x, y, z coordinates and an additional array to fill.

	This function is compliant with OGSTM-BFM code.
*/
vtkFloatArray *createVTKvecf3(const char *name, int nx, int ny, int nz, 
	double *a1, double *a2, double *a3) {

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
					// Project velocity on the staggered mesh
					double aux_u = (ii == 0) ? a1[ARRIND(ii,jj,kk,nx,ny)] : 
						(a1[ARRIND(ii-1,jj,kk,nx,ny)] + a1[ARRIND(ii,jj,kk,nx,ny)])/2.;
					double aux_v = (jj == 0) ? a2[ARRIND(ii,jj,kk,nx,ny)] : 
						(a2[ARRIND(ii,jj-1,kk,nx,ny)] + a2[ARRIND(ii,jj-1,kk,nx,ny)])/2.;
					double aux_w = (kk == 0) ? a3[ARRIND(ii,jj,kk,nx,ny)] : 
						(a3[ARRIND(ii,jj,kk-1,nx,ny)] + a3[ARRIND(ii,jj,kk,nx,ny)])/2.;
					// Set array value
					vtkArray->SetTuple3(VTKIND(ii,jj,kk,nx,ny),
						aux_u, aux_v, aux_w);
				}
			}
		}
	}

	return vtkArray;
}

/* CREATEVTKVECF3
	
	Creates a vtkFloatArray for a 3 component vector field given the array name, 
	its dimensions and an additional array to fill.

	This function is not compliant with OGSTM-BFM code.
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

/* GETCELLCOORDINATES

	Computes the cell centers given a VTKDataSet
*/
vtkFloatArray *getCellCoordinates(const char *varname, vtkDataSet *mesh) {

	int ncells = mesh->GetNumberOfCells();

	// Create a vector array to store the cell centers
	vtkFloatArray *vtkArray = createVTKvecf3(varname,ncells,NULL,NULL,NULL);

	// Loop the mesh and get the cell coordinates
	for (int cellId = 0; cellId < ncells; cellId++) {
		// Get number of points per cell
		int ncellp = mesh->GetCell(cellId)->GetNumberOfPoints();
		// Loop the points in the cell
		double xyzcen[3] = {0.,0.,0.};
		for (int pId = 0; pId < ncellp; pId++) {
			double xyz[3];
			mesh->GetCell(cellId)->GetPoints()->GetPoint(pId,xyz);
			// Average
			for (int ii = 0; ii < 3; ii++) xyzcen[ii] += xyz[ii]/ncellp;
		}
		// Set the Array
		vtkArray->SetTuple(cellId,xyzcen);
	}
	return vtkArray;
}

/* GETPOINTCOORDINATES

	Computes the point centers given a VTKDataSet
*/
vtkFloatArray *getPointCoordinates(const char *varname, vtkDataSet *mesh) {

	int npoints = mesh->GetNumberOfPoints();

	// Create a vector array to store the cell centers
	vtkFloatArray *vtkArray = createVTKvecf3(varname,npoints,NULL,NULL,NULL);

	// Loop the mesh and get the cell coordinates
	for (int pId = 0; pId < npoints; pId++) {
		// Get point
		double xyz[3];
		mesh->GetPoint(pId,xyz);
		// Set the Array
		vtkArray->SetTuple(pId,xyz);
	}
	return vtkArray;
}

/* COUNTUNIQUEZ
	
	Counts the number of unique values in Z direction
	and returns the values.
*/
int countUniqueZ(vtkFloatArray *vtkArray, int n, double epsi, 
	int *uniqueid, std::vector<double> &uniquevals) {

    int count = 0;
    uniquevals.clear();

    for (int ii = 0; ii < n; ii++) {
        double xyz[3];
        vtkArray->GetTuple(ii,xyz);

        int isunique = 1;
        for (int jj = 0; isunique && jj < ii; jj++) {
        	double xyz2[3];
        	vtkArray->GetTuple(jj,xyz2);
        	if ( fabs(xyz[2] - xyz2[2]) < epsi ) isunique = 0;
        }

        if ( isunique ) {
        	count++; uniquevals.push_back(xyz[2]);
        }
        uniqueid[ii] = -1;
        for (int jj = 0; jj < count; jj++) {
        	if ( fabs(xyz[2] - uniquevals.at(jj)) < epsi ) uniqueid[ii] = jj;
        }
    }
    
    return count;
}

/* COUNTUNIQUECOORDS
	
	Counts the number of unique coordinates in z direction and returns
	the true number of unique coordinates.
*/
int countUniqueCoords(vtkFloatArray *vtkArray, int n, int nunique, double epsi, 
	 int *uniqueid, std::vector<double> &uniquevals) {

	int count  = 0;

    for (int ii = 0; ii < n; ii++) {
    	// Recover a point in the mesh
        double xyz[3];
        vtkArray->GetTuple(ii,xyz);

        // Loop on the unique values to set the uniqueid
        uniqueid[ii] = -1;
        for (int jj = 0; jj < nunique-1 && uniqueid[ii] < 0; jj++) {
        	if ( uniquevals.at(jj) < xyz[2] )
        		uniqueid[ii] = jj;
        }
        // If the point hasn't been set, default it to the last one
        if (uniqueid[ii] < 0) uniqueid[ii] = nunique - 1;
        // Set count 
        count = (uniqueid[ii] > count-1) ? uniqueid[ii]+1 : count;
    }

    return count;
}

/* VTKGRAD3XY2

	Second order, face centered, approximation of the gradient for
	a vtk vectorial field.

	Output deri is in the form: [dudx, dudy, dudz, 
	                             dvdx, dvdy, dvdz, 
	                             dwdx, dwdy, dwdz]
*/
void vtkGrad3XY2(int ii, int jj, int kk, int nx, int ny, int nz,
	vtkFloatArray *vtkvecf, vtkFloatArray *vtkpoints, double deri[9]) {

	int ind = 0, ind1 = 0;
	double xyz[3], xyz1[3], val[3], val1[3];

	/*
		DERIVATIVES WITH RESPECT TO X
	*/
	if (ii == 0) {
		ind  = VTKIND(ii+1,jj,kk,nx,ny); ind1 = VTKIND(ii,jj,kk,nx,ny);
		ind  = VTKIND(ii+1,jj,kk,nx,ny); ind1 = VTKIND(ii,jj,kk,nx,ny);
	} else if (ii == nx-1) {
		ind  = VTKIND(ii,jj,kk,nx,ny); ind1 = VTKIND(ii-1,jj,kk,nx,ny);
		ind  = VTKIND(ii,jj,kk,nx,ny); ind1 = VTKIND(ii-1,jj,kk,nx,ny);
	} else {
		ind  = VTKIND(ii+1,jj,kk,nx,ny); ind1 = VTKIND(ii-1,jj,kk,nx,ny);
		ind  = VTKIND(ii+1,jj,kk,nx,ny); ind1 = VTKIND(ii-1,jj,kk,nx,ny);
	}
	// Recover values from vtkArrays
	vtkvecf  ->GetTuple(ind,val); vtkvecf  ->GetTuple(ind1,val1);
	vtkpoints->GetTuple(ind,xyz); vtkpoints->GetTuple(ind1,xyz1);
	// Compute the derivatives with respect to x
	deri[0] = (val[0] - val1[0])/(xyz[0] - xyz1[0]); // dudx
	deri[3] = (val[1] - val1[1])/(xyz[0] - xyz1[0]); // dvdx
	deri[6] = (val[2] - val1[2])/(xyz[0] - xyz1[0]); // dwdx

	/*
		DERIVATIVES WITH RESPECT TO Y
	*/
	if (jj == 0) {
		ind  = VTKIND(ii,jj+1,kk,nx,ny); ind1 = VTKIND(ii,jj,kk,nx,ny);
		ind  = VTKIND(ii,jj+1,kk,nx,ny); ind1 = VTKIND(ii,jj,kk,nx,ny);
	} else if (jj == ny-1) {
		ind  = VTKIND(ii,jj,kk,nx,ny); ind1 = VTKIND(ii,jj-1,kk,nx,ny);
		ind  = VTKIND(ii,jj,kk,nx,ny); ind1 = VTKIND(ii,jj-1,kk,nx,ny);
	} else {
		ind  = VTKIND(ii,jj+1,kk,nx,ny); ind1 = VTKIND(ii,jj-1,kk,nx,ny);
		ind  = VTKIND(ii,jj+1,kk,nx,ny); ind1 = VTKIND(ii,jj-1,kk,nx,ny);
	}
	// Recover values from vtkArrays
	vtkvecf  ->GetTuple(ind,val); vtkvecf  ->GetTuple(ind1,val1);
	vtkpoints->GetTuple(ind,xyz); vtkpoints->GetTuple(ind1,xyz1);
	// Compute the derivatives with respect to x
	deri[1] = (val[0] - val1[0])/(xyz[1] - xyz1[1]); // dudy
	deri[4] = (val[1] - val1[1])/(xyz[1] - xyz1[1]); // dvdy
	deri[7] = (val[2] - val1[2])/(xyz[1] - xyz1[1]); // dwdy

	/*
		DERIVATIVES WITH RESPECT TO Z
	*/
	if (kk == 0) {
		ind  = VTKIND(ii,jj,kk+1,nx,ny); ind1 = VTKIND(ii,jj,kk,nx,ny);
		ind  = VTKIND(ii,jj,kk+1,nx,ny); ind1 = VTKIND(ii,jj,kk,nx,ny);
	} else if (kk == nz-1) {
		ind  = VTKIND(ii,jj,kk,nx,ny); ind1 = VTKIND(ii,jj,kk-1,nx,ny);
		ind  = VTKIND(ii,jj,kk,nx,ny); ind1 = VTKIND(ii,jj,kk-1,nx,ny);
	} else {
		ind  = VTKIND(ii,jj,kk+1,nx,ny); ind1 = VTKIND(ii,jj,kk-1,nx,ny);
		ind  = VTKIND(ii,jj,kk+1,nx,ny); ind1 = VTKIND(ii,jj,kk-1,nx,ny);
	}
	// Recover values from vtkArrays
	vtkvecf  ->GetTuple(ind,val); vtkvecf  ->GetTuple(ind1,val1);
	vtkpoints->GetTuple(ind,xyz); vtkpoints->GetTuple(ind1,xyz1);
	// Compute the derivatives with respect to x
	deri[2] = (val[0] - val1[0])/(xyz[2] - xyz1[2]);  // dudz
	deri[5] = (val[1] - val1[1])/(xyz[2] - xyz1[2]);  // dvdz
	deri[8] = (val[2] - val1[2])/(xyz[2] - xyz1[2]);  // dwdz
}

/* VTKGRAD3XY4

	Fourth order, face centered, approximation of the gradient for
	a vtk vectorial field.

	Output deri is in the form: [dudx, dudy, dudz, 
	                             dvdx, dvdy, dvdz, 
	                             dwdx, dwdy, dwdz]
*/
void vtkGrad3XY4(int ii, int jj, int kk, int nx, int ny, int nz,
	vtkFloatArray *vtkvecf, vtkFloatArray *vtkpoints, double deri[9]) {
	
	double xyz[3], xyz1[3], val[3], val1[3], val2[3], val3[3];

	/*
		DERIVATIVES WITH RESPECT TO X
	*/
	if (ii <= 1) {
		vtkvecf->GetTuple(VTKIND(ii+2,jj,kk,nx,ny)  , val );
		vtkvecf->GetTuple(VTKIND(ii+1,jj,kk,nx,ny)  , val1);
		vtkvecf->GetTuple(VTKIND(ii,jj,kk,nx,ny)    , val2);
		vtkpoints->GetTuple(VTKIND(ii+1,jj,kk,nx,ny), xyz );
		vtkpoints->GetTuple(VTKIND(ii,jj,kk,nx,ny)  , xyz1);
		
		deri[0] = (-val[0] + 4*val1[0] - 3*val2[0])/2./(xyz[0] - xyz1[0]);  // dudx
		deri[3] = (-val[1] + 4*val1[1] - 3*val2[1])/2./(xyz[0] - xyz1[0]);  // dvdx
		deri[6] = (-val[2] + 4*val1[2] - 3*val2[2])/2./(xyz[0] - xyz1[0]);  // dwdx
	} else if (ii >= nx-2) {
		vtkvecf->GetTuple(VTKIND(ii,jj,kk,nx,ny)    , val );
		vtkvecf->GetTuple(VTKIND(ii-1,jj,kk,nx,ny)  , val1);
		vtkvecf->GetTuple(VTKIND(ii-2,jj,kk,nx,ny)  , val2);
		vtkpoints->GetTuple(VTKIND(ii,jj,kk,nx,ny)  , xyz );
		vtkpoints->GetTuple(VTKIND(ii-1,jj,kk,nx,ny), xyz1);

		deri[0] = (3*val[0] - 4*val1[0] + val2[0])/2./(xyz[0] - xyz1[0]); // dudx
		deri[3] = (3*val[1] - 4*val1[1] + val2[1])/2./(xyz[0] - xyz1[0]); // dvdx
		deri[6] = (3*val[2] - 4*val1[2] + val2[2])/2./(xyz[0] - xyz1[0]); // dwdx
	} else {
		vtkvecf->GetTuple(VTKIND(ii+2,jj,kk,nx,ny)  , val );
		vtkvecf->GetTuple(VTKIND(ii+1,jj,kk,nx,ny)  , val1);
		vtkvecf->GetTuple(VTKIND(ii-1,jj,kk,nx,ny)  , val2);
		vtkvecf->GetTuple(VTKIND(ii-2,jj,kk,nx,ny)  , val3);
		vtkpoints->GetTuple(VTKIND(ii+1,jj,kk,nx,ny), xyz );
		vtkpoints->GetTuple(VTKIND(ii-1,jj,kk,nx,ny), xyz1);

		deri[0] = (-val[0] + 8*val1[0] - 8*val2[0] + val3[0])/6./(xyz[0] - xyz1[0]); // dudx
		deri[3] = (-val[1] + 8*val1[1] - 8*val2[1] + val3[1])/6./(xyz[0] - xyz1[0]); // dvdx
		deri[6] = (-val[2] + 8*val1[2] - 8*val2[2] + val3[2])/6./(xyz[0] - xyz1[0]); // dwdx
	}
					
	/*
		DERIVATIVES WITH RESPECT TO Y
	*/
	if (jj <= 1) {
		vtkvecf->GetTuple(VTKIND(ii,jj+2,kk,nx,ny)  , val );
		vtkvecf->GetTuple(VTKIND(ii,jj+1,kk,nx,ny)  , val1);
		vtkvecf->GetTuple(VTKIND(ii,jj,kk,nx,ny)    , val2);
		vtkpoints->GetTuple(VTKIND(ii,jj+1,kk,nx,ny), xyz );
		vtkpoints->GetTuple(VTKIND(ii,jj,kk,nx,ny)  , xyz1);
		
		deri[1] = (-val[0] + 4*val1[0] - 3*val2[0])/2./(xyz[1] - xyz1[1]); // dudy
		deri[4] = (-val[1] + 4*val1[1] - 3*val2[1])/2./(xyz[1] - xyz1[1]); // dvdy
		deri[7] = (-val[2] + 4*val1[2] - 3*val2[2])/2./(xyz[1] - xyz1[1]); // dwdy
	} else if (jj >= ny-2) {
		vtkvecf->GetTuple(VTKIND(ii,jj,kk,nx,ny)    , val );
		vtkvecf->GetTuple(VTKIND(ii,jj-1,kk,nx,ny)  , val1);
		vtkvecf->GetTuple(VTKIND(ii,jj-2,kk,nx,ny)  , val2);
		vtkpoints->GetTuple(VTKIND(ii,jj,kk,nx,ny)  , xyz );
		vtkpoints->GetTuple(VTKIND(ii,jj-1,kk,nx,ny), xyz1);

		deri[1] = (3*val[0] - 4*val1[0] + val2[0])/2./(xyz[1] - xyz1[1]); // dudy
		deri[4] = (3*val[1] - 4*val1[1] + val2[1])/2./(xyz[1] - xyz1[1]); // dvdy
		deri[7] = (3*val[2] - 4*val1[2] + val2[2])/2./(xyz[1] - xyz1[1]); // dwdy
	} else {
		vtkvecf->GetTuple(VTKIND(ii,jj+2,kk,nx,ny)  , val );
		vtkvecf->GetTuple(VTKIND(ii,jj+1,kk,nx,ny)  , val1);
		vtkvecf->GetTuple(VTKIND(ii,jj-1,kk,nx,ny)  , val2);
		vtkvecf->GetTuple(VTKIND(ii,jj-2,kk,nx,ny)  , val3);
		vtkpoints->GetTuple(VTKIND(ii,jj+1,kk,nx,ny), xyz );
		vtkpoints->GetTuple(VTKIND(ii,jj-1,kk,nx,ny), xyz1);

		deri[1] = (-val[0] + 8*val1[0] - 8*val2[0] + val3[0])/6./(xyz[1] - xyz1[1]); // dudy
		deri[4] = (-val[1] + 8*val1[1] - 8*val2[1] + val3[1])/6./(xyz[1] - xyz1[1]); // dvdy
		deri[7] = (-val[2] + 8*val1[2] - 8*val2[2] + val3[2])/6./(xyz[1] - xyz1[1]); // dwdy
	}
	
	/*
		DERIVATIVES WITH RESPECT TO Z
	*/
	if (kk <= 1) {
		vtkvecf->GetTuple(VTKIND(ii,jj,kk+2,nx,ny)  , val );
		vtkvecf->GetTuple(VTKIND(ii,jj,kk+1,nx,ny)  , val1);
		vtkvecf->GetTuple(VTKIND(ii,jj,kk,nx,ny)    , val2);
		vtkpoints->GetTuple(VTKIND(ii,jj,kk+1,nx,ny), xyz );
		vtkpoints->GetTuple(VTKIND(ii,jj,kk,nx,ny)  , xyz1);
		
		deri[2] = (-val[0] + 4*val1[0] - 3*val2[0])/2./(xyz[2] - xyz1[2]); // dudz
		deri[5] = (-val[1] + 4*val1[1] - 3*val2[1])/2./(xyz[2] - xyz1[2]); // dvdz
		deri[8] = (-val[2] + 4*val1[2] - 3*val2[2])/2./(xyz[2] - xyz1[2]); // dwdz
	} else if (kk >= nz-2) {
		vtkvecf->GetTuple(VTKIND(ii,jj,kk,nx,ny)    , val );
		vtkvecf->GetTuple(VTKIND(ii,jj,kk-1,nx,ny)  , val1);
		vtkvecf->GetTuple(VTKIND(ii,jj,kk-2,nx,ny)  , val2);
		vtkpoints->GetTuple(VTKIND(ii,jj,kk,nx,ny)  , xyz );
		vtkpoints->GetTuple(VTKIND(ii,jj,kk-1,nx,ny), xyz1);

		deri[1] = (3*val[0] - 4*val1[0] + val2[0])/2./(xyz[2] - xyz1[2]); // dudz
		deri[4] = (3*val[1] - 4*val1[1] + val2[1])/2./(xyz[2] - xyz1[2]); // dvdz
		deri[7] = (3*val[2] - 4*val1[2] + val2[2])/2./(xyz[2] - xyz1[2]); // dwdz
	} else {
		vtkvecf->GetTuple(VTKIND(ii,jj,kk+2,nx,ny)  , val );
		vtkvecf->GetTuple(VTKIND(ii,jj,kk+1,nx,ny)  , val1);
		vtkvecf->GetTuple(VTKIND(ii,jj,kk-1,nx,ny)  , val2);
		vtkvecf->GetTuple(VTKIND(ii,jj,kk-2,nx,ny)  , val3);
		vtkpoints->GetTuple(VTKIND(ii,jj,kk+1,nx,ny), xyz );
		vtkpoints->GetTuple(VTKIND(ii,jj,kk-1,nx,ny), xyz1);

		deri[1] = (-val[0] + 8*val1[0] - 8*val2[0] + val3[0])/6./(xyz[2] - xyz1[2]); // dudz
		deri[4] = (-val[1] + 8*val1[1] - 8*val2[1] + val3[1])/6./(xyz[2] - xyz1[2]); // dvdz
		deri[7] = (-val[2] + 8*val1[2] - 8*val2[2] + val3[2])/6./(xyz[2] - xyz1[2]); // dwdz
	}
}

