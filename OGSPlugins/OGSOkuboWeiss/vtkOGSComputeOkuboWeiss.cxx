/*=========================================================================

  Program:   OGSComputeOkuboWeiss
  Module:    vtkOGSComputeOkuboWeiss.cxx

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkFloatArray.h"
#include "vtkCellData.h"
#include "vtkDataArraySelection.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkRectilinearGrid.h"

#include "vtkOGSComputeOkuboWeiss.h"

#include "vtkObjectFactory.h"

#define INDEX(ii,jj,kk,nx,ny) ( (nx-1)*(ny-1)*(kk) + (nx-1)*(jj) + (ii) )
#define GETVTKVAL1(vtkarray,ii,jj,kk,nx,ny)     ( (vtkarray)->GetTuple1((nx-1)*(ny-1)*(kk) + (nx-1)*(jj) + (ii)) )
#define GETVTKVAL3(vtkarray,ii,jj,kk,nx,ny)     ( (vtkarray)->GetTuple3((nx-1)*(ny-1)*(kk) + (nx-1)*(jj) + (ii)) )
#define SETVTKVAL1(vtkarray,val,ii,jj,kk,nx,ny) ( (vtkarray)->SetTuple1((nx-1)*(ny-1)*(kk) + (nx-1)*(jj) + (ii),(val)) )

vtkStandardNewMacro(vtkOGSComputeOkuboWeiss);

//----------------------------------------------------------------------------
void vtk3Grad2XY(int ii, int jj, int kk, int nx, int ny,
	vtkFloatArray *vtkvecf, vtkFloatArray *vtkx, vtkFloatArray *vtky, double deri[4]) {
/*
	Second order approximation for the gradient of a vtk vectorial
	field only in the x and y coordinates.

	Input deri is in the form: [dudx, dudy, dvdx, dvdy]
*/
	double val[3], val1[3], cor, cor1;

	// Bounds for X
	if (ii == 0) {
		vtkvecf->GetTuple(INDEX(ii+1,jj,kk,nx,ny),val);
		vtkvecf->GetTuple(INDEX(ii,jj,kk,nx,ny),val1);
		cor  = vtkx->GetTuple1(ii+1);
		cor1 = vtkx->GetTuple1(ii);  
	} else if (ii == nx-1) {
		vtkvecf->GetTuple(INDEX(ii,jj,kk,nx,ny),val);
		vtkvecf->GetTuple(INDEX(ii-1,jj,kk,nx,ny),val1);
		cor  = vtkx->GetTuple1(ii);
		cor1 = vtkx->GetTuple1(ii-1);
	} else {
		vtkvecf->GetTuple(INDEX(ii+1,jj,kk,nx,ny),val);
		vtkvecf->GetTuple(INDEX(ii-1,jj,kk,nx,ny),val1);
		cor  = vtkx->GetTuple1(ii+1);
		cor1 = vtkx->GetTuple1(ii-1);
	}
	deri[0] = (val[0] - val1[0])/(cor - cor1);
	deri[2] = (val[1] - val1[1])/(cor - cor1);
					
	// Bounds for y
	if (jj == 0) {
		vtkvecf->GetTuple(INDEX(ii,jj+1,kk,nx,ny),val);
		vtkvecf->GetTuple(INDEX(ii,jj,kk,nx,ny),val1);
		cor  = vtky->GetTuple1(jj+1);
		cor1 = vtky->GetTuple1(jj);
	} else if (jj == ny-1) {
		vtkvecf->GetTuple(INDEX(ii,jj,kk,nx,ny),val);
		vtkvecf->GetTuple(INDEX(ii,jj-1,kk,nx,ny),val1);
		cor  = vtky->GetTuple1(jj);
		cor1 = vtky->GetTuple1(jj-1);				
	} else {
		vtkvecf->GetTuple(INDEX(ii,jj+1,kk,nx,ny),val);
		vtkvecf->GetTuple(INDEX(ii,jj-1,kk,nx,ny),val1);
		cor  = vtky->GetTuple1(jj+1);
		cor1 = vtky->GetTuple1(jj-1);
	}
	deri[1] = (val[0] - val1[0])/(cor - cor1);
	deri[4] = (val[1] - val1[1])/(cor - cor1);
}

void vtk3Grad4XY(int ii, int jj, int kk, int nx, int ny,
	vtkFloatArray *vtkvecf, vtkFloatArray *vtkx, vtkFloatArray *vtky, double deri[4]) {
/*
	Fourth order approximation for the gradient of a vtk vectorial
	field only in the x and y coordinates.

	Input deri is in the form: [dudx, dudy, dvdx, dvdy]
*/
	double val[3], val1[3], val2[3], val3[3], cor, cor1;

	// Bounds for X
	if (ii <= 1) {
		vtkvecf->GetTuple(INDEX(ii+2,jj,kk,nx,ny),val);
		vtkvecf->GetTuple(INDEX(ii+1,jj,kk,nx,ny),val1);
		vtkvecf->GetTuple(INDEX(ii,jj,kk,nx,ny),val2);
		cor  = vtkx->GetTuple1(ii+1);
		cor1 = vtkx->GetTuple1(ii);

		deri[0] = (-val[0] + 4*val1[0] - 3*val2[0])/2/(cor - cor1);
		deri[2] = (-val[1] + 4*val1[1] - 3*val2[1])/2/(cor - cor1);	
	} else if (ii >= nx-2) {
		vtkvecf->GetTuple(INDEX(ii,jj,kk,nx,ny),val);
		vtkvecf->GetTuple(INDEX(ii-1,jj,kk,nx,ny),val1);
		vtkvecf->GetTuple(INDEX(ii-2,jj,kk,nx,ny),val2);
		cor  = vtkx->GetTuple1(ii);
		cor1 = vtkx->GetTuple1(ii-1);

		deri[0] = (3*val[0] - 4*val1[0] + val2[0])/2/(cor - cor1);
		deri[2] = (3*val[1] - 4*val1[1] + val2[1])/2/(cor - cor1);
	} else {
		vtkvecf->GetTuple(INDEX(ii+2,jj,kk,nx,ny),val);
		vtkvecf->GetTuple(INDEX(ii+1,jj,kk,nx,ny),val1);
		vtkvecf->GetTuple(INDEX(ii-1,jj,kk,nx,ny),val2);
		vtkvecf->GetTuple(INDEX(ii-2,jj,kk,nx,ny),val3);
		cor  = vtkx->GetTuple1(ii+1);
		cor1 = vtkx->GetTuple1(ii-1);

		deri[0] = (-val[0] + 8*val1[0] - 8*val2[0] + val3[0])/6/(cor - cor1);
		deri[2] = (-val[1] + 8*val1[1] - 8*val2[1] + val3[1])/6/(cor - cor1);
	}
					
	// Bounds for y
	if (jj <= 1) {
		vtkvecf->GetTuple(INDEX(ii,jj+2,kk,nx,ny),val);
		vtkvecf->GetTuple(INDEX(ii,jj+1,kk,nx,ny),val1);
		vtkvecf->GetTuple(INDEX(ii,jj,kk,nx,ny),val2);
		cor  = vtkx->GetTuple1(ii+1);
		cor1 = vtkx->GetTuple1(ii);

		deri[1] = (-val[0] + 4*val1[0] - 3*val2[0])/2/(cor - cor1);
		deri[4] = (-val[1] + 4*val1[1] - 3*val2[1])/2/(cor - cor1);	
	} else if (jj >= ny-2) {
		vtkvecf->GetTuple(INDEX(ii,jj,kk,nx,ny),val);
		vtkvecf->GetTuple(INDEX(ii,jj-1,kk,nx,ny),val1);
		vtkvecf->GetTuple(INDEX(ii,jj-2,kk,nx,ny),val2);
		cor  = vtkx->GetTuple1(ii);
		cor1 = vtkx->GetTuple1(ii-1);

		deri[1] = (3*val[0] - 4*val1[0] + val2[0])/2/(cor - cor1);
		deri[4] = (3*val[1] - 4*val1[1] + val2[1])/2/(cor - cor1);			
	} else {
		vtkvecf->GetTuple(INDEX(ii,jj+2,kk,nx,ny),val);
		vtkvecf->GetTuple(INDEX(ii,jj+1,kk,nx,ny),val1);
		vtkvecf->GetTuple(INDEX(ii,jj-1,kk,nx,ny),val2);
		vtkvecf->GetTuple(INDEX(ii,jj-2,kk,nx,ny),val3);
		cor  = vtkx->GetTuple1(ii+1);
		cor1 = vtkx->GetTuple1(ii-1);

		deri[1] = (-val[0] + 8*val1[0] - 8*val2[0] + val3[0])/6/(cor - cor1);
		deri[4] = (-val[1] + 8*val1[1] - 8*val2[1] + val3[1])/6/(cor - cor1);
	}
}

//----------------------------------------------------------------------------
vtkOGSComputeOkuboWeiss::vtkOGSComputeOkuboWeiss()
{
	this->field     = NULL;
	this->coef      = 0.2;
	this->grad_type = 0;
}

//----------------------------------------------------------------------------
vtkOGSComputeOkuboWeiss::~vtkOGSComputeOkuboWeiss()
{
	this->Setfield(NULL);
}

//----------------------------------------------------------------------------
int vtkOGSComputeOkuboWeiss::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
	// Get the info objects
	vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
	vtkInformation *outInfo = outputVector->GetInformationObject(0);

	// Get the input and output
	vtkRectilinearGrid *input = vtkRectilinearGrid::SafeDownCast(
		inInfo->Get(vtkDataObject::DATA_OBJECT()));
	vtkRectilinearGrid *output = vtkRectilinearGrid::SafeDownCast(
		outInfo->Get(vtkDataObject::DATA_OBJECT()));

	this->UpdateProgress(0.);

	// Recover XYZ coordinates
	vtkFloatArray *vtkxcoord = vtkFloatArray::SafeDownCast(
		input->GetXCoordinates());
	int nx = vtkxcoord->GetNumberOfTuples();
	vtkFloatArray *vtkycoord = vtkFloatArray::SafeDownCast(
		input->GetYCoordinates());
	int ny = vtkycoord->GetNumberOfTuples();
	vtkFloatArray *vtkzcoord = vtkFloatArray::SafeDownCast(
		input->GetZCoordinates());
	int nz = vtkzcoord->GetNumberOfTuples();

	// Recover velocity vector
	vtkFloatArray *vtkVeloc = vtkFloatArray::SafeDownCast(
		input->GetCellData()->GetArray(this->field));

	// Generate a new vtkFloatArray to store the Okubo-Weiss
	vtkFloatArray *vtkOW = vtkFloatArray::New();
	vtkOW->SetName("Okubo-Weiss");
	vtkOW->SetNumberOfComponents(1);   // Scalar field
	vtkOW->SetNumberOfTuples(nx*ny*nz);

	vtkFloatArray *vtkOWm = vtkFloatArray::New();
	vtkOWm->SetName("Okubo-Weiss mask");
	vtkOWm->SetNumberOfComponents(1);   // Scalar field
	vtkOWm->SetNumberOfTuples(nx*ny*nz);

	this->UpdateProgress(0.1);

	// Loop the mesh
	double OW_mean = 0., OW_std = 0.;

	for (int kk = 0; kk < nz; kk++) {
		// Computation of the Okubo-Weiss criterion at the surface
		if (kk == 0) {
			// First loop the 2D mesh and compute the gradients and mean
			for (int jj = 0; jj < ny; jj++) {
				for (int ii = 0; ii < nx; ii++) {
					double deri[4];
					// Selection of the gradient method
					switch (this->grad_type) {
						case 0:
							vtk3Grad2XY(ii,jj,kk,nx,ny,vtkVeloc,vtkxcoord,vtkycoord,deri);
							break;
						case 1:
							vtk3Grad4XY(ii,jj,kk,nx,ny,vtkVeloc,vtkxcoord,vtkycoord,deri);
							break;
					}
					// Rates of strain and stress
					double Sn = (deri[1] - deri[2]);
					double Ss = (deri[4] + deri[0]);
					double W  = (deri[4] - deri[0]);
					// Compute Okubo-Weiss
					double OW = Sn*Sn + Ss*Ss - W*W;
					// Compute the mean
					OW_mean += OW/nx/ny;
					// Set Okubo-Weiss
					SETVTKVAL1(vtkOW,OW,ii,jj,kk,nx,ny);
				}
			}
			// Up to this point the Okubo-Weiss criterion in the surface is computed
			// and stored in vtkOW and the mean in OW_mean.
			//
			// Now we work out the standard deviation
			for(int jj = 0; jj < ny; jj++){
				for (int ii = 0; ii < nx; ii++) {
					double aux = (GETVTKVAL1(vtkOW,ii,jj,kk,nx,ny) - OW_mean);
					OW_std += aux*aux/nx/ny;
				}
			}
			OW_std = this->coef*sqrt(OW_std);

			// We now have the standard deviation, let us work out the mask for
			// the surface
			for(int jj = 0; jj < ny; jj++){
				for (int ii = 0; ii < nx; ii++) {
					// Preallocate to zero
					SETVTKVAL1(vtkOWm,0.,ii,jj,kk,nx,ny);
					// Vorticity-dominated flow
					if (GETVTKVAL1(vtkOW,ii,jj,kk,nx,ny) < -OW_std)
						SETVTKVAL1(vtkOWm,-1.,ii,jj,kk,nx,ny);
					// Strain-dominated flow
					if (GETVTKVAL1(vtkOW,ii,jj,kk,nx,ny) > OW_std)
						SETVTKVAL1(vtkOWm,1.,ii,jj,kk,nx,ny);
				}
			}
		// Here kk > 1, we copy the values from the surface
		} else {
			for(int jj = 0; jj < ny; jj++){
				for (int ii = 0; ii < nx; ii++) {
					SETVTKVAL1(vtkOW,GETVTKVAL1(vtkOW,ii,jj,0,nx,ny),ii,jj,kk,nx,ny);
					SETVTKVAL1(vtkOWm,GETVTKVAL1(vtkOWm,ii,jj,0,nx,ny),ii,jj,kk,nx,ny);
				}
			}
		}
		this->UpdateProgress(0.1+0.8/nz*kk);
	}

	// Add arrays to input
	input->GetCellData()->AddArray(vtkOW);  vtkOW->Delete();
	input->GetCellData()->AddArray(vtkOWm); vtkOWm->Delete();

	// Copy the input grid
	this->UpdateProgress(1.);
	output->ShallowCopy(input);

	return 1;
}