/*=========================================================================

  Program:   OGSSelectCoast
  Module:    vtkOGSSelectCoast.cxx

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkTypeUInt8Array.h"
#include "vtkFloatArray.h"
#include "vtkDoubleArray.h"
#include "vtkCellData.h"
#include "vtkPointData.h"
#include "vtkDataArraySelection.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkUnstructuredGrid.h"

#include "vtkOGSSelectCoast.h"

#include "vtkObjectFactory.h"

#include <cstdint>

vtkStandardNewMacro(vtkOGSSelectCoast);

//----------------------------------------------------------------------------

/*
	Macro to set the array precision 
*/
#define FLDARRAY uint8_t
#define VTKARRAY vtkTypeUInt8Array

#include "../_utils/field.h"
#include "../_utils/vtkFields.hpp"

//----------------------------------------------------------------------------
void addCoasts(vtkDataArraySelection *CoastsDataArraySelection) {
	CoastsDataArraySelection->AddArray("Continental shelf");
	CoastsDataArraySelection->AddArray("Open Sea");
}

//----------------------------------------------------------------------------
vtkOGSSelectCoast::vtkOGSSelectCoast() {

	// Add the sub basins into the array
	this->CoastsDataArraySelection = vtkDataArraySelection::New();
	addCoasts(this->CoastsDataArraySelection);

	this->mask_field = NULL;
}

//----------------------------------------------------------------------------
vtkOGSSelectCoast::~vtkOGSSelectCoast()
{
	this->CoastsDataArraySelection->Delete();
	this->Setmask_field(NULL);
}

//----------------------------------------------------------------------------
void vtkOGSSelectCoast::PrintSelf(ostream& os, vtkIndent indent) {
	this->Superclass::PrintSelf(os, indent);
}

//----------------------------------------------------------------------------
int vtkOGSSelectCoast::RequestData(
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

	// Get the basins mask
	VTKARRAY *vtkmask = NULL;
	bool iscelld = true;

	vtkmask = VTKARRAY::SafeDownCast(input->GetCellData()->GetArray(this->mask_field));

	if (vtkmask == NULL) {
		vtkmask = VTKARRAY::SafeDownCast(input->GetPointData()->GetArray(this->mask_field));
		iscelld = false;
	}

	// Recover basins mask as a field
	field::Field<FLDARRAY> mask = VTK::createFieldfromVTK<VTKARRAY,FLDARRAY>(vtkmask);

	// Generate a new field (initialized at zero) that will be used as cutting mask
	field::Field<FLDARRAY> cutmask(mask.get_n(),1);

	this->UpdateProgress(0.2);

	// Loop and update cutting mask (Mesh loop, can be parallelized)
	field::Field<FLDARRAY>::iterator it_mask, it_cut;
	for (it_mask = mask.begin(), it_cut = cutmask.begin(); it_mask != mask.end(); ++it_mask,++it_cut) {
		it_cut[0] = 0;
		// Loop on the basins array selection
		for (int bid=0; bid < this->GetNumberOfCoastsArrays(); ++bid)
			// Test to zero is allowed when working with integers
			if ( this->GetCoastsArrayStatus(this->GetCoastsArrayName(bid)) && (it_mask[0] - (bid+1)) == 0 ) {
				it_cut[0] = 1; continue;
			}
	}

	// Convert field to vtkArray and add it to input
	VTKARRAY *vtkcutmask;
	vtkcutmask = VTK::createVTKfromField<VTKARRAY,FLDARRAY>("CutMask",cutmask);

	if (iscelld)
		input->GetCellData()->AddArray(vtkcutmask);
	else
		input->GetPointData()->AddArray(vtkcutmask);

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
	if (iscelld) {
		output->GetCellData()->RemoveArray("CutMask");
		output->GetCellData()->RemoveArray(this->mask_field);
	} else {
		output->GetPointData()->RemoveArray("CutMask");
		output->GetPointData()->RemoveArray(this->mask_field);		
	}
	vtkcutmask->Delete();

	// Return
	this->UpdateProgress(1.0);
	return 1;
}


//----------------------------------------------------------------------------
void vtkOGSSelectCoast::DisableAllCoastsArrays()
{
	this->CoastsDataArraySelection->DisableAllArrays();
}

void vtkOGSSelectCoast::EnableAllCoastsArrays()
{
	this->CoastsDataArraySelection->EnableAllArrays();
}

int vtkOGSSelectCoast::GetNumberOfCoastsArrays()
{
	return this->CoastsDataArraySelection->GetNumberOfArrays();
}

const char* vtkOGSSelectCoast::GetCoastsArrayName(int index)
{
	if (index >= (int)this->GetNumberOfCoastsArrays() || index < 0)
		return NULL;
	else
		return this->CoastsDataArraySelection->GetArrayName(index);
}

int vtkOGSSelectCoast::GetCoastsArrayIndex(const char* name)
{
	return this->CoastsDataArraySelection->GetArrayIndex(name);
}

int vtkOGSSelectCoast::GetCoastsArrayStatus(const char* name)
{
	return this->CoastsDataArraySelection->ArrayIsEnabled(name);
}

void vtkOGSSelectCoast::SetCoastsArrayStatus(const char* name, int status)
{
	if (status)
		this->CoastsDataArraySelection->EnableArray(name);
	else
		this->CoastsDataArraySelection->DisableArray(name);

	this->Modified();
}