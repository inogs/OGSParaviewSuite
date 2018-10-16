/*=========================================================================

  Program:   OGSComputeQCriterion
  Module:    vtkOGSComputeQCriterion.cxx

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

#include "vtkOGSComputeQCriterion.h"

#include "vtkObjectFactory.h"

namespace VTK
{
// Include the VTK Operations
#include "../_utils/vtkOperations.cpp"
}

vtkStandardNewMacro(vtkOGSComputeQCriterion);

//----------------------------------------------------------------------------
vtkOGSComputeQCriterion::vtkOGSComputeQCriterion()
{
	this->field     = NULL;
	this->coef      = 0.2;
	this->grad_type = 0;
}

//----------------------------------------------------------------------------
vtkOGSComputeQCriterion::~vtkOGSComputeQCriterion()
{
	this->Setfield(NULL);
}

//----------------------------------------------------------------------------
int vtkOGSComputeQCriterion::RequestData(
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

	output->ShallowCopy(input);

	this->UpdateProgress(0.);

	// Recover XYZ number of nodes
	int nx = input->GetXCoordinates()->GetNumberOfTuples();
	int ny = input->GetYCoordinates()->GetNumberOfTuples();
	int nz = input->GetZCoordinates()->GetNumberOfTuples();

	// Recover a mask
	// This is useful to perform a true average and srd on the Okubo-Weiss
	// Either the basins or the coast mask will work
	vtkFloatArray *vtkMask = vtkFloatArray::SafeDownCast(
		input->GetCellData()->GetArray("coast mask"));
	if (!vtkMask)
		vtkMask = vtkFloatArray::SafeDownCast(
		input->GetCellData()->GetArray("basins mask"));
	if (!vtkMask)
		vtkWarningMacro("Cannot find a working mask!\nResults might be affected.");

	// Recover weights
	vtkFloatArray *vtke1 = vtkFloatArray::SafeDownCast(
		input->GetCellData()->GetArray("e1"));
	vtkFloatArray *vtke2 = vtkFloatArray::SafeDownCast(
		input->GetCellData()->GetArray("e2"));
	vtkFloatArray *vtke3 = vtkFloatArray::SafeDownCast(
		input->GetCellData()->GetArray("e3"));

	if (this->grad_type > 1 && (vtke1 == NULL || vtke2 == NULL || vtke3 == NULL))
		vtkErrorMacro("Mesh weights (e1, e2 and e3) need to be loaded to proceed!");

	// Recover velocity vector
	vtkFloatArray *vtkVeloc = vtkFloatArray::SafeDownCast(
		input->GetCellData()->GetArray(this->field));

	// Compute the cell centers
	vtkFloatArray *vtkCellCenters = VTK::getCellCoordinates("Cell Centers",input);

	// Generate a new scalar array to store the Okubo-Weiss and the mask
	vtkFloatArray *vtkQ  = VTK::createVTKscaf("Q-criterion",nx-1,ny-1,nz-1,NULL);
	vtkFloatArray *vtkQm = VTK::createVTKscaf("Q-criterion mask",nx-1,ny-1,nz-1,NULL);

	this->UpdateProgress(0.1);

	// Loop the mesh
	double Q_mean = 0., sum_weights = 0., Q_std = 0.;
	double deri[9]     = {0.,0.,0.,0.,0.,0.,0.,0.,0.};
	double deri_old[9] = {0.,0.,0.,0.,0.,0.,0.,0.,0.};

	// Loop on the 3D mesh and compute Q-criterion
	for (int kk = 0; kk < nz-1; kk++) {
		for (int jj = 0; jj < ny-1; jj++) {
			for (int ii = 0; ii < nx-1; ii++) {
				// Check the mask, skip the point if we are out of the mask
				if (vtkMask && vtkMask->GetTuple1(CLLIND(ii,jj,kk,nx,ny)) < 0)
					continue;
				// Store the previous derivative
				if (this->grad_type > 1) {
					for (int dd = 0; dd < 9; dd++) 
						deri_old[dd] = deri[dd];
				}
				// Recover e1u, e2v and e3w
				double e1[4], e2[4], e3[4];
				vtke1->GetTuple(CLLIND(ii,jj,kk,nx,ny),e1);
				vtke2->GetTuple(CLLIND(ii,jj,kk,nx,ny),e2);
				vtke3->GetTuple(CLLIND(ii,jj,kk,nx,ny),e3);
				// Selection of the gradient method
				switch (this->grad_type) {
					case 0: // Second order, face centered gradient
							// This gradient is unsafe as it relies on the mesh projection
						VTK::vtkGrad3XY2(ii,jj,kk,nx-1,ny-1,nz-1,vtkVeloc,vtkCellCenters,deri);
						break;
					case 1: // Fourth order, face centered gradient
							// This gradient is unsafe as it relies on the mesh projection
						VTK::vtkGrad3XY4(ii,jj,kk,nx-1,ny-1,nz-1,vtkVeloc,vtkCellCenters,deri);
						break;
					case 2:	// OGSTM-BFM approach according to the NEMO handbook
							// This gradient is safe as it relies on the code implementation
						VTK::vtkGrad3OGS1(ii,jj,kk,nx-1,ny-1,nz-1,
							vtkVeloc,e1[1],e2[2],e3[3],deri);
						break;
					case 3:	// 2nd order OGSTM-BFM approach
							// This gradient is experimental
						VTK::vtkGrad3OGS2(ii,jj,kk,nx-1,ny-1,nz-1,
							vtkVeloc,e1[1],e2[2],e3[3],deri);
						break;
					case 4:	// 4th order OGSTM-BFM approach
							// This gradient is experimental
						VTK::vtkGrad3OGS4(ii,jj,kk,nx-1,ny-1,nz-1,
							vtkVeloc,e1[1],e2[2],e3[3],deri);
						break;
				}
				// Interpolate on the centered grid
				if (this->grad_type > 1) {
					// Interpolation first component
					VTK::facecen2cellcen(ii,jj,kk,nx-1,ny-1,
						deri_old,deri,vtke1,vtke2,vtke3,deri);
					// Interpolation second component
					VTK::facecen2cellcen(ii,jj,kk,nx-1,ny-1,
						deri_old+3,deri+3,vtke1,vtke2,vtke3,deri+3);
					// Interpolation third component
					VTK::facecen2cellcen(ii,jj,kk,nx-1,ny-1,
						deri_old+6,deri+6,vtke1,vtke2,vtke3,deri+6);		
				}
				// Compute Q-criterion
				double Q = -0.5*(deri[0]*deri[0] + deri[4]*deri[4] + deri[8]*deri[8])
						-deri[1]*deri[3]-deri[2]*deri[6]-deri[5]*deri[7];
				Q *= -2;
				// Compute the mean
				Q_mean      += Q*e1[0]*e2[0];
				sum_weights += e1[0]*e2[0];
				// Set Okubo-Weiss
				vtkQ->SetTuple1(CLLIND(ii,jj,kk,nx,ny),Q);
			}
		}
		this->UpdateProgress(0.1+0.3/(nz-1.)*kk);
	}
	// Up to this point the Okubo-Weiss criterion in the surface is computed
	// and stored in vtkOW and the mean in OW_mean.
	Q_mean /= sum_weights;

	// Now we work out the standard deviation
	for (int kk = 0; kk < nz-1; kk++) {
		for(int jj = 0; jj < ny-1; jj++){
			for (int ii = 0; ii < nx-1; ii++) {
				// Check the mask, skip the point if we are out of the mask
				if (vtkMask && vtkMask->GetTuple1(CLLIND(ii,jj,kk,nx,ny)) < 0)
					continue;
				double e1[4], e2[4];
				vtke1->GetTuple(CLLIND(ii,jj,kk,nx,ny),e1);
				vtke2->GetTuple(CLLIND(ii,jj,kk,nx,ny),e2);
				double Q  = vtkQ->GetTuple1(CLLIND(ii,jj,kk,nx,ny));
				Q_std    += (Q - Q_mean)*(Q - Q_mean)*e1[0]*e2[0];
			}
		}
		this->UpdateProgress(0.4+0.3/(nz-1.)*kk);
	}
	Q_std = this->coef*sqrt(Q_std/sum_weights);

	// We now have the standard deviation, let us work out the mask for
	// the surface
	for (int kk = 0; kk < nz-1; kk++) {
		for(int jj = 0; jj < ny-1; jj++){
			for (int ii = 0; ii < nx-1; ii++) {
				// Check the mask, skip the point if we are out of the mask
				if (vtkMask && vtkMask->GetTuple1(CLLIND(ii,jj,kk,nx,ny)) < 0)
					continue;
				double Q  = vtkQ->GetTuple1(CLLIND(ii,jj,kk,nx,ny));
				// Vorticity-dominated flow
				if (Q < -Q_std)
					vtkQm->SetTuple1(CLLIND(ii,jj,kk,nx,ny),-1.);
				// Strain-dominated flow
				if (Q > Q_std)
					vtkQm->SetTuple1(CLLIND(ii,jj,kk,nx,ny),1.);
			}
		}
		this->UpdateProgress(0.7+0.3/(nz-1.)*kk);
	}
		
	vtkCellCenters->Delete();

	// Add arrays to input
	output->GetCellData()->AddArray(vtkQ);  vtkQ->Delete();
	output->GetCellData()->AddArray(vtkQm); vtkQm->Delete();

	// Copy the input grid
	this->UpdateProgress(1.);
	return 1;
}