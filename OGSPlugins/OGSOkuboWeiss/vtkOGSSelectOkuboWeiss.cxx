/*=========================================================================

  Program:   OGSSelectBasin
  Module:    vtkOGSSelectOkuboWeiss.cxx

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkFloatArray.h"
#include "vtkCellData.h"
#include "vtkPointData.h"
#include "vtkDataArraySelection.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkUnstructuredGrid.h"

#include "vtkOGSSelectOkuboWeiss.h"

#include "vtkObjectFactory.h"

vtkStandardNewMacro(vtkOGSSelectOkuboWeiss);

//----------------------------------------------------------------------------
vtkOGSSelectOkuboWeiss::vtkOGSSelectOkuboWeiss() {

	// Add the sub basins into the array
	this->OWDataArraySelection = vtkDataArraySelection::New();
	this->OWDataArraySelection->AddArray("Vorticity dominated");
	this->OWDataArraySelection->AddArray("Strain dominated");
	this->OWDataArraySelection->AddArray("Background field");

	this->mask_field = NULL;
}

//----------------------------------------------------------------------------
vtkOGSSelectOkuboWeiss::~vtkOGSSelectOkuboWeiss()
{
	this->OWDataArraySelection->Delete();
	this->Setmask_field(NULL);
}

//----------------------------------------------------------------------------
void vtkOGSSelectOkuboWeiss::PrintSelf(ostream& os, vtkIndent indent) {
	this->Superclass::PrintSelf(os, indent);
}

//----------------------------------------------------------------------------
int vtkOGSSelectOkuboWeiss::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
	// Get the info objects
	vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
	vtkInformation *outInfo = outputVector->GetInformationObject(0);

	// Get the input and output
	vtkDataSet *input = vtkDataSet::SafeDownCast(
		inInfo->Get(vtkDataObject::DATA_OBJECT()));
	vtkUnstructuredGrid *output = vtkUnstructuredGrid::SafeDownCast(
		outInfo->Get(vtkDataObject::DATA_OBJECT()));

	this->UpdateProgress(0.0);

	// Assume we are working with cell arrays
	vtkFloatArray *ow_mask = vtkFloatArray::SafeDownCast(
		input->GetCellData()->GetArray(this->mask_field));
	int npoints = input->GetNumberOfCells();

	// Check if we are really working with cell arrays
	if (ow_mask == NULL) {
		vtkFloatArray *ow_mask = vtkFloatArray::SafeDownCast(
			input->GetPointData()->GetArray(this->mask_field));
		int npoints = input->GetNumberOfPoints();	
	}	

	// Generate a new vtkFloatArray that will be used for the mask
	vtkFloatArray *mask = vtkFloatArray::New();
	mask->SetName("CutMask");
	mask->SetNumberOfComponents(1);   // Scalar field
	mask->SetNumberOfTuples(npoints);

	this->UpdateProgress(0.2);

	// Loop the mesh and set the mask
	for (int ii = 0; ii < npoints; ii++) {
		// Recover value from basins mask
		int owmask_val = (int)(ow_mask->GetTuple1(ii));
		// Initialize mask to zero
		mask->SetTuple1(ii,0);
		// Set the conditions
		if (this->GetOWArrayStatus("Vorticity dominated") && owmask_val == -1)
				mask->SetTuple1(ii,1);
		if (this->GetOWArrayStatus("Strain dominated") && owmask_val == 1)
				mask->SetTuple1(ii,1);
		if (this->GetOWArrayStatus("Background field") && owmask_val == 0)
				mask->SetTuple1(ii,1);
	}

	// Add array to input
	input->GetCellData()->AddArray(mask);

	this->UpdateProgress(0.4);

	// Force ThresholdBetween to obtain values that are greater than 0
	this->Superclass::ThresholdBetween(0.5,1.);

	// Force to use the CutMask to produce the Threshold
	this->Superclass::SetInputArrayToProcess(0,0,0,
		vtkDataObject::FIELD_ASSOCIATION_CELLS,"CutMask");

	this->UpdateProgress(0.6);

	// Run the actual threshold filter
	this->Superclass::RequestData(NULL,inputVector,outputVector);

	this->UpdateProgress(0.8);

	// Cleanup the output by deleting the CutMask and the basins mask
	output->GetCellData()->RemoveArray("CutMask");
	output->GetCellData()->RemoveArray(this->mask_field);

	// Cleanup
	mask->Delete(); 

	// Return
	this->UpdateProgress(1.0);
	return 1;
}

//----------------------------------------------------------------------------
void vtkOGSSelectOkuboWeiss::DisableAllOWArrays()
{
	this->OWDataArraySelection->DisableAllArrays();
}

void vtkOGSSelectOkuboWeiss::EnableAllOWArrays()
{
	this->OWDataArraySelection->EnableAllArrays();
}

int vtkOGSSelectOkuboWeiss::GetNumberOfOWArrays()
{
	return this->OWDataArraySelection->GetNumberOfArrays();
}

const char* vtkOGSSelectOkuboWeiss::GetOWArrayName(int index)
{
	if (index >= (int)this->GetNumberOfOWArrays() || index < 0)
		return NULL;
	else
		return this->OWDataArraySelection->GetArrayName(index);
}

int vtkOGSSelectOkuboWeiss::GetOWArrayIndex(const char* name)
{
	return this->OWDataArraySelection->GetArrayIndex(name);
}

int vtkOGSSelectOkuboWeiss::GetOWArrayStatus(const char* name)
{
	return this->OWDataArraySelection->ArrayIsEnabled(name);
}

void vtkOGSSelectOkuboWeiss::SetOWArrayStatus(const char* name, int status)
{
	if (status)
		this->OWDataArraySelection->EnableArray(name);
	else
		this->OWDataArraySelection->DisableArray(name);

	this->Modified();
}