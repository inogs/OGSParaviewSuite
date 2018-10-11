/*=========================================================================

  Program:   OGSComputeOkuboWeiss
  Module:    vtkOGSComputeOkuboWeiss.cxx

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
#include "vtkCellData.h"
#include "vtkDataArraySelection.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkRectilinearGrid.h"

#include "vtkOGSComputeOkuboWeiss.h"

#include "vtkObjectFactory.h"

namespace VTK
{
// Include the VTK Operations
#include "../_utils/vtkOperations.cpp"
}

//#define INDEX(ii,jj,kk,nx,ny) ( (nx-1)*(ny-1)*(kk) + (nx-1)*(jj) + (ii) )
//#define GETVTKVAL1(vtkarray,ii,jj,kk,nx,ny)     ( (vtkarray)->GetTuple1((nx-1)*(ny-1)*(kk) + (nx-1)*(jj) + (ii)) )
//#define GETVTKVAL3(vtkarray,ii,jj,kk,nx,ny)     ( (vtkarray)->GetTuple3((nx-1)*(ny-1)*(kk) + (nx-1)*(jj) + (ii)) )
//#define SETVTKVAL1(vtkarray,val,ii,jj,kk,nx,ny) ( (vtkarray)->SetTuple1((nx-1)*(ny-1)*(kk) + (nx-1)*(jj) + (ii),(val)) )

vtkStandardNewMacro(vtkOGSComputeOkuboWeiss);

//----------------------------------------------------------------------------




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

	// Recover XYZ number of nodes
	int nx = input->GetXCoordinates()->GetNumberOfTuples();
	int ny = input->GetXCoordinates()->GetNumberOfTuples();
	int nz = input->GetXCoordinates()->GetNumberOfTuples();

	// Recover weights
	vtkFloatArray *vtke1t = vtkFloatArray::SafeDownCast(
		input->GetCellData()->GetArray("e1t"));
	vtkFloatArray *vtke2t = vtkFloatArray::SafeDownCast(
		input->GetCellData()->GetArray("e2t"));
	vtkFloatArray *vtke1u = vtkFloatArray::SafeDownCast(
		input->GetCellData()->GetArray("e1u"));
	vtkFloatArray *vtke2v = vtkFloatArray::SafeDownCast(
		input->GetCellData()->GetArray("e2v"));
	vtkFloatArray *vtke3w = vtkFloatArray::SafeDownCast(
		input->GetCellData()->GetArray("e3w"));

	// Recover velocity vector
	vtkFloatArray *vtkVeloc = vtkFloatArray::SafeDownCast(
		input->GetCellData()->GetArray(this->field));

	// Compute the cell centers
	vtkFloatArray *vtkCellCenters = VTK::getCellCoordinates("Cell Centers",input);

	// Generate a new scalar array to store the Okubo-Weiss and the mask
	vtkFloatArray *vtkOW  = VTK::createVTKscaf("Okubo-Weiss",nx*ny*nz,NULL);
	vtkFloatArray *vtkOWm = VTK::createVTKscaf("Okubo-Weiss mask",nx*ny*nz,NULL);

	this->UpdateProgress(0.1);

	// Loop the mesh
	double OW_mean = 0., sum_weights = 0., OW_std = 0.;
	double deri_old[9] = {0.,0.,0.,0.,0.,0.,0.,0.,0.};

	for (int kk = 0; kk < nz; kk++) {
		// Computation of the Okubo-Weiss criterion at the surface
		if (kk == 0) {
			// First loop the 2D mesh and compute the gradients and mean
			for (int jj = 0; jj < ny; jj++) {
				for (int ii = 0; ii < nx; ii++) {
					double deri[9] = {0.,0.,0.,0.,0.,0.,0.,0.,0.}; // We will not use the z derivatives
					// Selection of the gradient method
					switch (this->grad_type) {
						case 0: // Second order, face centered gradient
							    // This gradient is unsafe as it relies on the mesh projection
							VTK::vtkGrad3XY2(ii,jj,kk,nx,ny,nz,vtkVeloc,vtkCellCenters,deri);
							break;
						case 1: // Fourth order, face centered gradient
							    // This gradient is unsafe as it relies on the mesh projection
							VTK::vtkGrad3XY4(ii,jj,kk,nx,ny,nz,vtkVeloc,vtkCellCenters,deri);
							break;
						case 2:	// OGSTM-BFM approach according to the NEMO handbook
								// This gradient is safe as it relies on the code implementation
							for (int dd = 0; dd < 9; dd++) deri_old[dd] = deri[dd];
							VTK::vtkGrad3OGS1(ii,jj,kk,nx,ny,nz,
								vtkVeloc,vtke1u,vtke2v,vtke3w,deri,deri_old);
							break;
					}
					// Rates of strain and stress
					double Sn = (deri[1] - deri[3]); // dudy - dvdx
					double Ss = (deri[4] + deri[0]); // dvdy + dudx
					double W  = (deri[4] - deri[0]); // dvdy - dudx
					// Compute Okubo-Weiss
					double OW = Sn*Sn + Ss*Ss - W*W;
					// Compute the mean
					double e1t   = vtke1t->GetTuple1(VTKIND(ii,jj,kk,nx,ny));
					double e2t   = vtke2t->GetTuple1(VTKIND(ii,jj,kk,nx,ny));
					OW_mean     += OW*e1t*e2t;
					sum_weights += e1t*e2t;
					// Set Okubo-Weiss
					vtkOW->SetTuple1(VTKIND(ii,jj,kk,nx,ny),OW);
				}
			}
			// Up to this point the Okubo-Weiss criterion in the surface is computed
			// and stored in vtkOW and the mean in OW_mean.
			OW_mean /= sum_weights;

			// Now we work out the standard deviation
			for(int jj = 0; jj < ny; jj++){
				for (int ii = 0; ii < nx; ii++) {
					double e1t = vtke1t->GetTuple1(VTKIND(ii,jj,kk,nx,ny));
					double e2t = vtke2t->GetTuple1(VTKIND(ii,jj,kk,nx,ny));
					double OW  = vtkOW->GetTuple1(VTKIND(ii,jj,kk,nx,ny));
					OW_std    += (OW - OW_mean)*(OW - OW_mean)*e1t*e2t;
				}
			}
			OW_std = this->coef*sqrt(OW_std/sum_weights);

			// We now have the standard deviation, let us work out the mask for
			// the surface
			for(int jj = 0; jj < ny; jj++){
				for (int ii = 0; ii < nx; ii++) {
					double OW  = vtkOW->GetTuple1(VTKIND(ii,jj,kk,nx,ny));
					// Vorticity-dominated flow
					if (OW < -OW_std)
						vtkOWm->SetTuple1(VTKIND(ii,jj,kk,nx,ny),-1.);
					// Strain-dominated flow
					if (OW > OW_std)
						vtkOWm->SetTuple1(VTKIND(ii,jj,kk,nx,ny),1.);
				}
			}
		// Here kk > 1, we copy the values from the surface
		} else {
			for(int jj = 0; jj < ny; jj++){
				for (int ii = 0; ii < nx; ii++) {
					vtkOW->SetTuple1(VTKIND(ii,jj,kk,nx,ny),
						vtkOW->GetTuple1(VTKIND(ii,jj,0,nx,ny)));
					vtkOWm->SetTuple1(VTKIND(ii,jj,kk,nx,ny),
						vtkOWm->GetTuple1(VTKIND(ii,jj,0,nx,ny)));
				}
			}
		}
		this->UpdateProgress(0.1+0.8/nz*kk);
	}
	vtkCellCenters->Delete();

	// Add arrays to input
	input->GetCellData()->AddArray(vtkOW);  vtkOW->Delete();
	input->GetCellData()->AddArray(vtkOWm); vtkOWm->Delete();

	// Copy the input grid
	this->UpdateProgress(1.);
	output->ShallowCopy(input);

	return 1;
}