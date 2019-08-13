/*=========================================================================

	Program:   OGSWriter
	Module:    vtkOGSNPZWriter.cxx

	Copyright (c) 2019 Arnau Miro, OGS
	All rights reserved.

		 This software is distributed WITHOUT ANY WARRANTY; without even
		 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
		 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkOGSNPZWriter.h"

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
#include <string>
#include <vtksys/SystemTools.hxx>

vtkStandardNewMacro(vtkOGSNPZWriter);

//----------------------------------------------------------------------------

/*
	Macro to set the array precision 
*/
#define FLDARRAY double
#define FLDMASK  uint8_t
#define VTKARRAY vtkDoubleArray
#define VTKMASK  vtkTypeUInt8Array

#include "../_utils/V3.h"
#include "../_utils/field.h"
#include "../_utils/vtkFields.hpp"
#include "../_utils/vtkOperations.hpp"

#include"cnpy.hpp"

//----------------------------------------------------------------------------
vtkOGSNPZWriter::vtkOGSNPZWriter() : FileName(nullptr), varname(nullptr), dfact(1000.), append(0), singlevar(0), ii_cur(0) {}

//----------------------------------------------------------------------------
vtkOGSNPZWriter::~vtkOGSNPZWriter() {
	this->SetFileName(nullptr);
	this->Setvarname(nullptr);
}

//----------------------------------------------------------------------------
int vtkOGSNPZWriter::FillInputPortInformation( int vtkNotUsed(port), vtkInformation* info) {
	info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
	return 1;
}

//----------------------------------------------------------------------------
int vtkOGSNPZWriter::Write() {
	// Make sure we only export from one mesh
	if(this->GetNumberOfInputConnections(0) != 1) {
		vtkErrorMacro("Exactly one input required."); 
		return 0;
	}

	// Always write even if the data hasn't changed
	this->Modified();
	this->Update();

	return 1;
}

//----------------------------------------------------------------------------
int vtkOGSNPZWriter::ProcessRequest(vtkInformation* request, 
	vtkInformationVector** inputVector, vtkInformationVector* outputVector) {
	return this->Superclass::ProcessRequest(request, inputVector, outputVector);
}

//----------------------------------------------------------------------------
int vtkOGSNPZWriter::RequestData(vtkInformation* request, 
	vtkInformationVector** inputVector, vtkInformationVector* vtkNotUsed(outputVector)) {

	// Tell the pipeline to start looping.
	if (this->ii_cur == this->ii_start && this->timeseries)
		request->Set(vtkStreamingDemandDrivenPipeline::CONTINUE_EXECUTING(), 1);

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
    	
    	// Format the filename
    	std::string path       = vtksys::SystemTools::GetFilenamePath(this->FileName);
    	std::string fnamenoext = vtksys::SystemTools::GetFilenameWithoutLastExtension(this->FileName);
    	std::string ext        = vtksys::SystemTools::GetFilenameLastExtension(this->FileName);

    	// File name
    	std::string fname = path + std::string("/") + fnamenoext + std::string(".") + std::string(buff) + ext;
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
void vtkOGSNPZWriter::WriteData() {
	// Recover the input
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

	// If we are not appending data, compute the points or the cell centers
	// and create the field to save the variables
	if (!this->append) {
		// Compute the mesh depending if we're dealing with point or cell data
		v3::V3v xyz = (iscelld) ? VTK::getVTKCellCenters(input,this->dfact) : VTK::getVTKCellPoints(input,this->dfact);
		// Store on the npz file
		size_t len = (size_t)(xyz.len());
		double *aux = new double[3*len]; aux = xyz.todouble();
		cnpy::npz_save<double>(this->FileName,"xyz",aux,{len,3},"w");
		delete [] aux;
	}

	// Save the variables, depending if we are saving a single variable or multiple variables
	if (this->singlevar) {
		// Load and store the variable
		if (std::string(this->varname) == std::string("basins mask") || 
			std::string(this->varname) == std::string("coast mask")) {
			// Use mask uint8
			field::Field<FLDMASK> array;
			array = VTK::createFieldfromVTK<VTKMASK,FLDMASK>( VTKMASK::SafeDownCast(vtkArray) );
			size_t Nn = (size_t)(array.get_n()), Nm = (size_t)(array.get_m());
			cnpy::npz_save<FLDMASK>(this->FileName,this->varname,array.data(),{Nn,Nm},"a");
		} else {
			// Use normal array
			field::Field<FLDARRAY> array;
			array = VTK::createFieldfromVTK<VTKARRAY,FLDARRAY>( VTKARRAY::SafeDownCast(vtkArray) );
			cnpy::npz_save<FLDARRAY>(this->FileName,this->varname,array.data(),{(size_t)(array.get_n()),(size_t)(array.get_m())},"a");
		}
	} else {
		int nvars = (iscelld) ? input->GetCellData()->GetNumberOfArrays() : input->GetPointData()->GetNumberOfArrays();
		for (int varId = 0; varId < nvars; ++varId) {
			// Load vtkArray
			vtkArray = (iscelld) ? input->GetCellData()->GetArray(varId) : input->GetPointData()->GetArray(varId);
			// Load and store the variable
			if (std::string(vtkArray->GetName()) == std::string("basins mask") || 
				std::string(vtkArray->GetName()) == std::string("coast mask")) {
				// Use mask uint8
				field::Field<FLDMASK> array;
				array = VTK::createFieldfromVTK<VTKMASK,FLDMASK>( VTKMASK::SafeDownCast(vtkArray) );
				cnpy::npz_save<FLDMASK>(this->FileName,vtkArray->GetName(),array.data(),{(size_t)(array.get_n()),(size_t)(array.get_m())},"a");
			} else {
				// Use normal array
				field::Field<FLDARRAY> array;
				array = VTK::createFieldfromVTK<VTKARRAY,FLDARRAY>( VTKARRAY::SafeDownCast(vtkArray) );
				cnpy::npz_save<FLDARRAY>(this->FileName,vtkArray->GetName(),array.data(),{(size_t)(array.get_n()),(size_t)(array.get_m())},"a");
			}
		}
	}
}

//----------------------------------------------------------------------------
void vtkOGSNPZWriter::SetStartEnd(const int val1, const int val2) {
	this->ii_start = val1;
	this->ii_end   = val2;
	this->ii_cur   = val1;
	this->Modified();
}

