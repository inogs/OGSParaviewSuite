/*=========================================================================

  Program:   OGSAnnotateDateTime
  Module:    vtkOGSAnnotateDateTime.cxx

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkOGSAnnotateDateTime.h"

#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#include "vtkDataSet.h"
#include "vtkFieldData.h"
#include "vtkStringArray.h"

#include "vtkObjectFactory.h"

#include <ctime>
#include <chrono>

vtkStandardNewMacro(vtkOGSAnnotateDateTime);

#ifdef PARAVIEW_USE_MPI
#include "vtkMultiProcessController.h"
vtkCxxSetObjectMacro(vtkOGSAnnotateDateTime, Controller, vtkMultiProcessController);
#endif

//----------------------------------------------------------------------------
vtkOGSAnnotateDateTime::vtkOGSAnnotateDateTime() {
	this->TimeFormat  = NULL;
	this->useMetadata = 0;
	this->nProcs      = 0;
	this->procId      = 0;

	#ifdef PARAVIEW_USE_MPI
		this->Controller = NULL;
		this->SetController(vtkMultiProcessController::GetGlobalController());
	#endif
}

//----------------------------------------------------------------------------
vtkOGSAnnotateDateTime::~vtkOGSAnnotateDateTime() {
	this->SetTimeFormat(0);

	#ifdef PARAVIEW_USE_MPI
		this->SetController(NULL);	
	#endif
}

//----------------------------------------------------------------------------
int vtkOGSAnnotateDateTime::RequestData(vtkInformation* request,
	vtkInformationVector** inputVector, vtkInformationVector* outputVector) {

	// Recover the Date string
	vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);

	// Stop all threads except from the master to execute
	#ifdef PARAVIEW_USE_MPI
	if (this->procId > 0) return 1;
	#endif

	// Decide how to write the time, whether to use metadata or to use a conversion
	// from the timestep to a string format
	struct tm tm = {0};
	char buff[256];

	if (this->useMetadata) {
		vtkDataSet *input = vtkDataSet::SafeDownCast(
			inInfo->Get(vtkDataObject::DATA_OBJECT()));

		// Recover metadata
		vtkStringArray *vtkmetadata = vtkStringArray::SafeDownCast(
			input->GetFieldData()->GetAbstractArray("Metadata"));

		// If successful, modify the entry by parsing the string
		// using the time library
		if (vtkmetadata) {
			strptime(vtkmetadata->GetValue(0).c_str(),"%Y%m%d-%H:%M:%S",&tm);
			strftime(buff,256,this->TimeFormat,&tm);
			this->Superclass::SetFormat(buff);
		}
	} else {
		// Recover the current timestep
		double requestedTimeValue = inInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());

		// Convert to struct tm
		time_t time = std::chrono::system_clock::to_time_t(std::chrono::system_clock::time_point(
			std::chrono::duration_cast<std::chrono::seconds>(std::chrono::duration<double>(requestedTimeValue))));
		tm = *localtime(&time);

		// Format and display
		strftime(buff,256,this->TimeFormat,&tm);
		this->Superclass::SetFormat(buff);
	}

	// Run RequestData from PythonAnnotationFilter
	return this->Superclass::RequestData(request,inputVector,outputVector);
}
