/*=========================================================================

	Program:   OGSWriter
	Module:    vtkOGSFieldWriter.cxx

	Copyright (c) 2019 Arnau Miro, OGS
	All rights reserved.

		 This software is distributed WITHOUT ANY WARRANTY; without even
		 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
		 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkOGSFieldWriter.h"

#include "vtkCell.h"
#include "vtkPoints.h"
#include "vtkCallbackCommand.h"
#include "vtkCommand.h"
#include "vtkExecutive.h"
#include "vtkCellData.h"
#include "vtkPointData.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkTypeUInt8Array.h"
#include "vtkFloatArray.h"
#include "vtkAbstractArray.h"
#include "vtkDoubleArray.h"
#include "vtkStringArray.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#include "vtkObjectFactory.h"

#include <cmath>
#include <ctime>
#include <chrono>
#include <vector>
#include <vtksys/SystemTools.hxx>

vtkStandardNewMacro(vtkOGSFieldWriter);

//----------------------------------------------------------------------------
#include "macros.h"
#include "V3.h"
#include "field.h"
#include "fieldOperations.h"
#include "vtkFields.h"
#include "vtkOperations.h"

//----------------------------------------------------------------------------
vtkOGSFieldWriter::vtkOGSFieldWriter() : FileName(nullptr), varname(nullptr), dfact(1000.), fstride(1) {}

//----------------------------------------------------------------------------
vtkOGSFieldWriter::~vtkOGSFieldWriter() {
	this->SetFileName(nullptr);
	this->Setvarname(nullptr);
}

//----------------------------------------------------------------------------
int vtkOGSFieldWriter::FillInputPortInformation( int vtkNotUsed(port), vtkInformation* info) {
	info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
	return 1;
}

//----------------------------------------------------------------------------
int vtkOGSFieldWriter::Write() {
	// Make sure we only export from one mesh
	if(this->GetNumberOfInputConnections(0) != 1) {
		vtkErrorMacro("Exactly one input required."); 
		return 0;
	}

	// Extract filename format
	this->path       = vtksys::SystemTools::GetFilenamePath(this->FileName);
	this->fnamenoext = vtksys::SystemTools::GetFilenameWithoutLastExtension(this->FileName);
	this->ext        = vtksys::SystemTools::GetFilenameLastExtension(this->FileName);

	// Always write even if the data hasn't changed
	this->Modified();
	this->Update();

	return 1;
}

//----------------------------------------------------------------------------
int vtkOGSFieldWriter::ProcessRequest(vtkInformation* request, 
	vtkInformationVector** inputVector, vtkInformationVector* outputVector) {
	return this->Superclass::ProcessRequest(request, inputVector, outputVector);
}

//----------------------------------------------------------------------------
int vtkOGSFieldWriter::RequestData(vtkInformation* request, 
	vtkInformationVector** inputVector, vtkInformationVector* vtkNotUsed(outputVector)) {

	// Tell the pipeline to start looping.
	if (this->ii_cur == this->ii_start && this->timeseries) {
		request->Set(vtkStreamingDemandDrivenPipeline::CONTINUE_EXECUTING(), 1);
		// Number of timesteps failsafe
		int ntsteps = inputVector[0]->GetInformationObject(0)->Length( vtkStreamingDemandDrivenPipeline::TIME_STEPS() );
		this->ii_end = (this->ii_end > ntsteps) ? ntsteps : this->ii_end;
	}

	// Handle the timestep
	struct tm tm = {0}; char buff[256];
	double *inTimes = inputVector[0]->GetInformationObject(0)->Get( vtkStreamingDemandDrivenPipeline::TIME_STEPS() );

	if (inTimes && this->timeseries) {
		// Recover the current timestep and set it
		double timeReq = inTimes[this->ii_cur];
		inputVector[0]->GetInformationObject(0)->Set( vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP(),timeReq );
		
		// Convert to struct tm
		time_t time = std::chrono::system_clock::to_time_t(std::chrono::system_clock::time_point(
			std::chrono::duration_cast<std::chrono::seconds>(std::chrono::duration<double>(timeReq))));
		tm = *localtime(&time);

		// Format the time
		strftime(buff,256,"%Y%m%d-%H:%M:%S",&tm);

		// File name
		std::string fname = this->path + std::string("/") + this->fnamenoext + std::string(".") + std::string(buff) + this->ext;
		this->SetFileName(fname.c_str());
	}	
	
	// Write the data
	this->WriteData();

	if (this->timeseries) {
		this->ii_cur++;
		if (this->ii_cur >= this->ii_end) {
			// Tell the pipeline to stop looping.
			request->Remove(vtkStreamingDemandDrivenPipeline::CONTINUE_EXECUTING());
			this->ii_cur = this->ii_start;
		}
	}

	return 1;
}

//----------------------------------------------------------------------------
void vtkOGSFieldWriter::WriteData() {
	// Make sure we only export from one mesh
	if(this->GetNumberOfInputConnections(0) != 1)
		vtkErrorMacro("Exactly one input required.");

	// Recover the input
	vtkInformation *inInfo = vtkInformation::SafeDownCast( this->GetExecutive()->GetInputInformation(0,0) );
	vtkDataSet *input = vtkDataSet::SafeDownCast( this->GetExecutive()->GetInputData(0,0) );

	// Check that varname is not metadata
	if (std::string(this->varname) == std::string("Metadata")) {
		vtkWarningMacro("Metadata cannot be the input variable! Aborting.");
		return;
	}

	// Recover Metadata array (depth factor)
	vtkStringArray *vtkmetadata = vtkStringArray::SafeDownCast(
		input->GetFieldData()->GetAbstractArray("Metadata"));
	this->dfact = (vtkmetadata != NULL) ? std::stod( vtkmetadata->GetValue(2) ) : this->dfact;
		
	if (vtkmetadata == NULL) 
		vtkWarningMacro("Field array Metadata not found! Using user input value.");

	// Try to load as Cell data
	// Use the variable in varname to check if we are using CellData or PointData
	bool iscelld = true; vtkAbstractArray *vtkArray;
	vtkArray = input->GetCellData()->GetArray(this->varname);
	// Then array might be point data
	if (!vtkArray) {
		vtkArray = input->GetPointData()->GetArray(this->varname);
		iscelld = false;
	}

	// Recover current timestep
	double tstep = inInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());

	// Compute the points or the cell centers
	v3::V3v xyz = (iscelld) ? VTK::getVTKCellCenters(input,this->dfact) : VTK::getVTKCellPoints(input,this->dfact);

	// Recover the field and write
	if (std::string(this->varname) == std::string("basins mask") || 
		std::string(this->varname) == std::string("coast mask")  ||
		std::string(this->varname) == std::string("land mask")) {
		// Use mask uint8
		field::Field<FLDMASK> array;
		array = VTK::createFieldfromVTK<VTKMASK,FLDMASK>( VTKMASK::SafeDownCast(vtkArray) );
		// Write field
		field::WriteField<FLDMASK>(this->FileName,(float)(tstep),xyz,array,(this->fstride) ? true : false);
	} else {
		// Use normal array
		field::Field<FLDARRAY> array;
		array = VTK::createFieldfromVTK<VTKARRAY,FLDARRAY>( VTKARRAY::SafeDownCast(vtkArray) );
		// Write field
		field::WriteField<FLDARRAY>(this->FileName,(float)(tstep),xyz,array,(this->fstride) ? true : false);
	}	
}

//----------------------------------------------------------------------------
void vtkOGSFieldWriter::SetStartEnd(const int val1, const int val2) {
	this->ii_start = val1;
	this->ii_end   = val2;
	this->ii_cur   = val1;
	this->Modified();
}