/*=========================================================================

  Program:   OGSComputeLambda2Criterion
  Module:    vtkOGSComputeLambda2Criterion.cxx

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

#include "vtkOGSComputeLambda2Criterion.h"

#include "vtkObjectFactory.h"

namespace VTK
{
// Include the VTK Operations
#include "../_utils/vtkOperations.cpp"
}
namespace MMATH
{
// Include matrix operations
#include "../_utils/matrixOperations.cpp"
}


vtkStandardNewMacro(vtkOGSComputeLambda2Criterion);

//----------------------------------------------------------------------------
vtkOGSComputeLambda2Criterion::vtkOGSComputeLambda2Criterion()
{
	this->field     = NULL;
	this->grad_type = 0;
}

//----------------------------------------------------------------------------
vtkOGSComputeLambda2Criterion::~vtkOGSComputeLambda2Criterion()
{
	this->Setfield(NULL);
}

//----------------------------------------------------------------------------
int vtkOGSComputeLambda2Criterion::RequestData(
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

	if (this->grad_type > 1 && (vtke1 == NULL || vtke2 == NULL || vtke3 == NULL)) {
		vtkErrorMacro("Mesh weights (e1, e2 and e3) need to be loaded to proceed!");
		return 0;
	}

	// Recover velocity vector
	vtkFloatArray *vtkVeloc = vtkFloatArray::SafeDownCast(
		input->GetCellData()->GetArray(this->field));

	// Compute the cell centers
	vtkFloatArray *vtkCellCenters = VTK::getCellCoordinates("Cell Centers",input);

	// Generate a new scalar array to store the Okubo-Weiss and the mask
	vtkFloatArray *vtkL2 = VTK::createVTKscaf("Lambda2 criterion",nx-1,ny-1,nz-1,NULL);

	this->UpdateProgress(0.1);

	// Loop the mesh
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

				// Compute symmetric tensor S
				double S[9];
				S[0] = 0.5*(deri[0] + deri[0]);
				S[1] = 0.5*(deri[1] + deri[3]);
				S[2] = 0.5*(deri[2] + deri[6]);
				S[3] = 0.5*(deri[3] + deri[1]);
				S[4] = 0.5*(deri[4] + deri[4]);
				S[5] = 0.5*(deri[5] + deri[7]);
				S[6] = 0.5*(deri[6] + deri[2]);
				S[7] = 0.5*(deri[7] + deri[5]);
				S[8] = 0.5*(deri[8] + deri[8]);
				
				// Compute asymetric tensor O
				double O[9];
				O[0] = 0.;
				O[1] = 0.5*(deri[1] - deri[3]);
				O[2] = 0.5*(deri[2] - deri[6]);
				O[3] = 0.5*(deri[3] - deri[1]);
				O[4] = 0.;
				O[5] = 0.5*(deri[5] - deri[7]);
				O[6] = 0.5*(deri[6] - deri[2]);
				O[7] = 0.5*(deri[7] - deri[5]);
				O[8] = 0.;
				
				// Compute R
				double R[9] = {0,0,0,0,0,0,0,0,0};

				MMATH::math_33_3_product_add(S, S, R);
				MMATH::math_33_3_product_add(O, O, R);
				
				// Find the eigenvalues of R
				double e[3] = {0,0,0}, v[9];
				int it_num, rot_num;

				MMATH::math_jacobi_eigenvalue(3, R, 100, v, e, it_num, rot_num);

				// Set Lambda2
				vtkL2->SetTuple1(CLLIND(ii,jj,kk,nx,ny),e[1]);
			}
		}
		this->UpdateProgress(0.1+0.9/(nz-1.)*kk);
	}
		
	vtkCellCenters->Delete();

	// Add arrays to input
	output->GetCellData()->AddArray(vtkL2);  vtkL2->Delete();

	// Copy the input grid
	this->UpdateProgress(1.);
	return 1;
}