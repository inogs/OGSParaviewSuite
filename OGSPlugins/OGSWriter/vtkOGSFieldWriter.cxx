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

#include<vector>
#include<string>

vtkStandardNewMacro(vtkOGSFieldWriter);

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
#include "../_utils/fieldOperations.hpp"
#include "../_utils/vtkFields.hpp"
#include "../_utils/vtkOperations.hpp"

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
		std::string(this->varname) == std::string("coast mask")) {
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
int vtkOGSFieldWriter::Write() {
	return Superclass::Write();
}