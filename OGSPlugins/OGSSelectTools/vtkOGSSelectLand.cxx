/*=========================================================================

  Program:   OGSSelectTools
  Module:    vtkOGSSelectLand.cxx

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkOGSSelectLand.h"

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
vtkCxxSetObjectMacro(vtkOGSSelectLand, Controller, vtkMultiProcessController);
#endif

vtkStandardNewMacro(vtkOGSSelectLand);

//----------------------------------------------------------------------------

#include "macros.h"
#include "field.h"
#include "vtkFields.h"

//----------------------------------------------------------------------------
vtkOGSSelectLand::vtkOGSSelectLand() {

	this->mask_field = NULL;
	this->nProcs     = 0;
	this->procId     = 0;

	#ifdef PARAVIEW_USE_MPI
		this->Controller = NULL;
		this->SetController(vtkMultiProcessController::GetGlobalController());
	#endif
}

//----------------------------------------------------------------------------
vtkOGSSelectLand::~vtkOGSSelectLand() {
	this->Setmask_field(NULL);

	#ifdef PARAVIEW_USE_MPI
		this->SetController(NULL);	
	#endif
}

//----------------------------------------------------------------------------
void vtkOGSSelectLand::PrintSelf(ostream& os, vtkIndent indent) {
	this->Superclass::PrintSelf(os, indent);
}

//----------------------------------------------------------------------------
int vtkOGSSelectLand::RequestData(vtkInformation *vtkNotUsed(request),
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

	// Decide whether we have cell or point data
	int n_cell_vars  = input->GetCellData()->GetNumberOfArrays();
	int n_point_vars = input->GetPointData()->GetNumberOfArrays();

	bool iscelld = (n_cell_vars > n_point_vars) ? true : false;

	if (iscelld) {
		// Force to use the CutMask to produce the Threshold
		this->Superclass::SetInputArrayToProcess(0,0,0,
			vtkDataObject::FIELD_ASSOCIATION_CELLS,this->mask_field);
	} else {
		// Force to use the CutMask to produce the Threshold
		this->Superclass::SetInputArrayToProcess(0,0,0,
			vtkDataObject::FIELD_ASSOCIATION_POINTS,this->mask_field);
	}

	this->UpdateProgress(0.4);

	// Force ThresholdBetween to obtain values that are greater than 0
	this->Superclass::ThresholdBetween(0.5,1.);

	this->UpdateProgress(0.6);

	// Run the actual threshold filter
	this->Superclass::RequestData(NULL,inputVector,outputVector);

	this->UpdateProgress(0.8);

	// Cleanup the output by deleting the CutMask and the basins mask
	if (iscelld) {
		output->GetCellData()->RemoveArray(this->mask_field);
	} else {
		output->GetPointData()->RemoveArray(this->mask_field);		
	}

	// Return
	this->UpdateProgress(1.0);
	return 1;
}
