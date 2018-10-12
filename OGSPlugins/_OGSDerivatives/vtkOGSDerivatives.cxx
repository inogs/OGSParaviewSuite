/*=========================================================================

  Program:   OGSDerivatives
  Module:    vtkOGSDerivatives.cxx

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

#include "vtkOGSDerivatives.h"

#include "vtkObjectFactory.h"

namespace VTK
{
// Include the VTK Operations
#include "../_utils/vtkOperations.cpp"
}

vtkStandardNewMacro(vtkOGSDerivatives);

//----------------------------------------------------------------------------
vtkOGSDerivatives::vtkOGSDerivatives()
{
	this->field     = NULL;
	this->grad_type = 0;
}

//----------------------------------------------------------------------------
vtkOGSDerivatives::~vtkOGSDerivatives()
{
	this->Setfield(NULL);
}

//----------------------------------------------------------------------------
int vtkOGSDerivatives::RequestData(
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

	// Try to load the array, if we fail it is a point array and we should
	// inform the user that point arrays are not treated in this function
	vtkFloatArray *vtkArray = vtkFloatArray::SafeDownCast(
		input->GetCellData()->GetArray(this->field));
	if (vtkArray == NULL)
		vtkErrorMacro("Input array "<<this->field<<"is not a cell array."<<
			"This function can only deal with cell arrays!");

	// Here we have a valid vtkArray to process. We must decide whether to act as
	// a scalar array or a vector array.
	int ncomp = vtkArray->GetNumberOfComponents();

	if (ncomp == 1)
		// Scalar array
		this->ScalarArrayDerivatives(vtkArray,output);
	else if (ncomp == 3)
		// Vector array
		this->VectorArrayDerivatives(vtkArray,output);
	else
		// I don't know what this is...
		vtkErrorMacro("Input array "<<this->field<<"is neither a scalar or vector array!");

	// Copy the input grid
	this->UpdateProgress(1.);
	return 1;
}

//----------------------------------------------------------------------------
void vtkOGSDerivatives::ScalarArrayDerivatives(vtkFloatArray *vtkArray, vtkRectilinearGrid *mesh) {
	// Recover XYZ number of nodes
	int nx = mesh->GetXCoordinates()->GetNumberOfTuples();
	int ny = mesh->GetYCoordinates()->GetNumberOfTuples();
	int nz = mesh->GetZCoordinates()->GetNumberOfTuples();

	// Recover weights
	vtkFloatArray *vtke1 = vtkFloatArray::SafeDownCast(
		mesh->GetCellData()->GetArray("e1"));
	vtkFloatArray *vtke2 = vtkFloatArray::SafeDownCast(
		mesh->GetCellData()->GetArray("e2"));
	vtkFloatArray *vtke3 = vtkFloatArray::SafeDownCast(
		mesh->GetCellData()->GetArray("e3"));

	if (vtke1 == NULL || vtke2 == NULL || vtke3 == NULL)
		vtkErrorMacro("Mesh weights (e1, e2 and e3) need to be loaded to proceed!");

	// Compute the cell centers
	vtkFloatArray *vtkCellCenters = VTK::getCellCoordinates("Cell Centers",mesh);

	// Preallocate output array, which will be a vector array
	char varname[256];
	sprintf(varname,"grad(%s)",vtkArray->GetName());
	vtkFloatArray *vtkGradArray = VTK::createVTKvecf3(varname,nx,ny,nz,NULL,NULL,NULL);

	// Loop the mesh cells
	double deri[3]     = {0.,0.,0.};
	double deri_old[3] = {0.,0.,0.};
	
	for (int kk = 0; kk < nz-1; kk++) {
		for (int jj = 0; jj < ny-1; jj++) {
			for (int ii = 0; ii < nx-1; ii++) {
				// Store the previous derivative
				if (this->grad_type > 1) {
					for (int dd = 0; dd < 3; dd++) 
						deri_old[dd] = deri[dd];
				}
				// Recover e1u, e2v and e3w
				double e1[4], e2[4], e3[4];
				vtke1->GetTuple(CLLIND(ii,jj,0,nx,ny),e1);
				vtke2->GetTuple(CLLIND(ii,jj,0,nx,ny),e2);
				vtke3->GetTuple(CLLIND(ii,jj,0,nx,ny),e3);
				// Selection of the gradient method
				switch (this->grad_type) {
					case 0: // Second order, face centered gradient
							// This gradient is unsafe as it relies on the mesh projection
						VTK::vtkGradXY2(ii,jj,kk,nx-1,ny-1,nz-1,vtkArray,vtkCellCenters,deri);
						break;
					case 1: // Fourth order, face centered gradient
							// This gradient is unsafe as it relies on the mesh projection
						VTK::vtkGradXY4(ii,jj,kk,nx-1,ny-1,nz-1,vtkArray,vtkCellCenters,deri);
						break;
					case 2:	// OGSTM-BFM approach according to the NEMO handbook
							// This gradient is safe as it relies on the code implementation
						VTK::vtkGradOGS1(ii,jj,kk,nx-1,ny-1,nz-1,
							vtkArray,e1[1],e2[2],e3[3],deri);
						break;
					case 3:	// 2nd order OGSTM-BFM approach
							// This gradient is experimental
						VTK::vtkGradOGS2(ii,jj,kk,nx-1,ny-1,nz-1,
							vtkArray,e1[1],e2[2],e3[3],deri);
						break;
					case 4:	// 4th order OGSTM-BFM approach
							// This gradient is experimental
						VTK::vtkGradOGS4(ii,jj,kk,nx-1,ny-1,nz-1,
							vtkArray,e1[1],e2[2],e3[3],deri);
						break;
				}
				// Interpolate on the centered grid
				if (this->grad_type > 1)
					VTK::facecen2cellcen(ii,jj,0,nx-1,ny-1,
						deri_old,deri,vtke1,vtke2,vtke3,deri);
				// Add to array
				vtkGradArray->SetTuple(CLLIND(ii,jj,kk,nx,ny),deri);
			}
		}
		this->UpdateProgress(0.+1./(nz-1.)*kk);
	}
	// Add array to mesh
	mesh->GetCellData()->AddArray(vtkGradArray); 
	vtkGradArray->Delete();
	vtkCellCenters->Delete();
}

//----------------------------------------------------------------------------
void vtkOGSDerivatives::VectorArrayDerivatives(vtkFloatArray *vtkArray, vtkRectilinearGrid *mesh) {


}