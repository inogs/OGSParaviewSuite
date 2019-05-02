/*=========================================================================

  Program:   OGSSelectBasin
  Module:    vtkOGSSelectBasin.cxx

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkOGSSelectBasin.h"

#include "vtkTypeUInt8Array.h"
#include "vtkFloatArray.h"
#include "vtkDoubleArray.h"
#include "vtkCellData.h"
#include "vtkPointData.h"
#include "vtkDataArraySelection.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkUnstructuredGrid.h"

#include "vtkObjectFactory.h"

#include <cstdint>

#ifdef __linux__
// Include OpenMP when working with GCC
#include <omp.h>
#define OMP_NUM_THREADS omp_get_num_threads()
#define OMP_THREAD_NUM  omp_get_thread_num()
#else
#define OMP_NUM_THREADS 1
#define OMP_THREAD_NUM  0
#endif

#ifdef PARAVIEW_USE_MPI
#include "vtkMultiProcessController.h"
vtkCxxSetObjectMacro(vtkOGSSelectBasin, Controller, vtkMultiProcessController);
#endif

vtkStandardNewMacro(vtkOGSSelectBasin);

//----------------------------------------------------------------------------

/*
	Macro to set the array precision 
*/
#define FLDMASK uint8_t
#define VTKMASK vtkTypeUInt8Array

#include "../_utils/field.h"
#include "../_utils/vtkFields.hpp"

//----------------------------------------------------------------------------
void addSubBasins(vtkDataArraySelection *BasinsDataArraySelection) {
	BasinsDataArraySelection->AddArray("Alboran Sea");
	BasinsDataArraySelection->AddArray("South Western Mediterranean (west)");
	BasinsDataArraySelection->AddArray("South Western Mediterranean (east)");
	BasinsDataArraySelection->AddArray("North Western Mediterranean");
	BasinsDataArraySelection->AddArray("Northern Tyrrhenian");
	BasinsDataArraySelection->AddArray("Southern Tyrrhenian");
	BasinsDataArraySelection->AddArray("Northern Adriatic");
	BasinsDataArraySelection->AddArray("Southern Adriatic");
	BasinsDataArraySelection->AddArray("Aegean Sea");
	BasinsDataArraySelection->AddArray("Western Ionian");
	BasinsDataArraySelection->AddArray("Eastern Ionian");
	BasinsDataArraySelection->AddArray("Northern Ionian");
	BasinsDataArraySelection->AddArray("Western Levantine");
	BasinsDataArraySelection->AddArray("Northern Levantine");
	BasinsDataArraySelection->AddArray("Southern Levantine");
	BasinsDataArraySelection->AddArray("Eastern Levantine");
}

//----------------------------------------------------------------------------
vtkOGSSelectBasin::vtkOGSSelectBasin() {
	// Add the sub basins into the array
	this->BasinsDataArraySelection = vtkDataArraySelection::New();
	addSubBasins(this->BasinsDataArraySelection);

	this->mask_field = NULL;
	this->nProcs     = 0;
	this->procId     = 0;

	#ifdef PARAVIEW_USE_MPI
		this->Controller = NULL;
		this->SetController(vtkMultiProcessController::GetGlobalController());
	#endif
}

//----------------------------------------------------------------------------
vtkOGSSelectBasin::~vtkOGSSelectBasin() {
	this->BasinsDataArraySelection->Delete();
	this->Setmask_field(NULL);

	#ifdef PARAVIEW_USE_MPI
		this->SetController(NULL);	
	#endif
}

//----------------------------------------------------------------------------
void vtkOGSSelectBasin::PrintSelf(ostream& os, vtkIndent indent) {
	this->Superclass::PrintSelf(os, indent);
}

//----------------------------------------------------------------------------
int vtkOGSSelectBasin::RequestData(vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector, vtkInformationVector *outputVector) {
	// Get the info objects
	vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
	vtkInformation *outInfo = outputVector->GetInformationObject(0);

	// Stop all threads except from the master to execute
	#ifdef PARAVIEW_USE_MPI
	if (this->procId > 0) return 1;
	#endif

	// Get the input and output
	vtkDataSet *input = vtkDataSet::SafeDownCast(
		inInfo->Get(vtkDataObject::DATA_OBJECT()));
	vtkUnstructuredGrid *output = vtkUnstructuredGrid::SafeDownCast(
		outInfo->Get(vtkDataObject::DATA_OBJECT()));

	this->UpdateProgress(0.0);

	// Get the basins mask
	VTKMASK *vtkmask = NULL;
	bool iscelld = true;

	vtkmask = VTKMASK::SafeDownCast(input->GetCellData()->GetArray(this->mask_field));

	if (vtkmask == NULL) {
		vtkmask = VTKMASK::SafeDownCast(input->GetPointData()->GetArray(this->mask_field));
		iscelld = false;
	}

	if (vtkmask == NULL) {
		vtkErrorMacro("Cannot load mask field "<<this->mask_field<<"!");
		return 0;
	}

	// Recover basins mask as a field
	field::Field<FLDMASK> mask = VTK::createFieldfromVTK<VTKMASK,FLDMASK>(vtkmask);

	// Generate a new field (initialized at zero) that will be used as cutting mask
	field::Field<FLDMASK> cutmask(mask.get_n(),1);

	this->UpdateProgress(0.2);

	// Loop and update cutting mask (Mesh loop, can be parallelized)
	#pragma omp parallel shared(mask,cutmask)
	{
	for (int ii = 0 + OMP_THREAD_NUM; ii < mask.get_n(); ii += OMP_NUM_THREADS) {
		cutmask[ii][0] = 0;
		// Loop on the basins array selection
		for (int bid=0; bid < this->GetNumberOfBasinsArrays(); ++bid)
			if (this->GetBasinsArrayStatus(this->GetBasinsArrayName(bid)) && mask[ii][bid]) { cutmask[ii][0] = 1; continue; }
	}
	}

	// Convert field to vtkArray and add it to input
	VTKMASK *vtkcutmask;
	vtkcutmask = VTK::createVTKfromField<VTKMASK,FLDMASK>("CutMask",cutmask);

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
void vtkOGSSelectBasin::DisableAllBasinsArrays()
{
	this->BasinsDataArraySelection->DisableAllArrays();
}

void vtkOGSSelectBasin::EnableAllBasinsArrays()
{
	this->BasinsDataArraySelection->EnableAllArrays();
}

int vtkOGSSelectBasin::GetNumberOfBasinsArrays()
{
	return this->BasinsDataArraySelection->GetNumberOfArrays();
}

const char* vtkOGSSelectBasin::GetBasinsArrayName(int index)
{
	if (index >= (int)this->GetNumberOfBasinsArrays() || index < 0)
		return NULL;
	else
		return this->BasinsDataArraySelection->GetArrayName(index);
}

int vtkOGSSelectBasin::GetBasinsArrayIndex(const char* name)
{
	return this->BasinsDataArraySelection->GetArrayIndex(name);
}

int vtkOGSSelectBasin::GetBasinsArrayStatus(const char* name)
{
	return this->BasinsDataArraySelection->ArrayIsEnabled(name);
}

void vtkOGSSelectBasin::SetBasinsArrayStatus(const char* name, int status)
{
	if (status)
		this->BasinsDataArraySelection->EnableArray(name);
	else
		this->BasinsDataArraySelection->DisableArray(name);

	this->Modified();
}