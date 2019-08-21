/*=========================================================================

  Program:   OGSDerivatives
  Module:    vtkOGSDerivatives.cxx

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

	 This software is distributed WITHOUT ANY WARRANTY; without even
	 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
	 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkOGSDerivatives.h"

#include "vtkCell.h"
#include "vtkPoints.h"
#include "vtkCellData.h"
#include "vtkPointData.h"
#include "vtkFieldData.h"
#include "vtkFloatArray.h"
#include "vtkDoubleArray.h"
#include "vtkCellData.h"
#include "vtkDataArraySelection.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkRectilinearGrid.h"

#include "vtkObjectFactory.h"

#include <string>

#ifdef PARAVIEW_USE_MPI
#include "vtkMultiProcessController.h"
vtkCxxSetObjectMacro(vtkOGSDerivatives, Controller, vtkMultiProcessController);
#endif

vtkStandardNewMacro(vtkOGSDerivatives);

//----------------------------------------------------------------------------

/*
	Macro to set the array precision 
*/
#define FLDARRAY double
#define VTKARRAY vtkDoubleArray

// V3.h and field.h defined in vtkOGSDerivatives.h
#include "../_utils/fieldOperations.hpp"
#include "../_utils/vtkFields.hpp"
#include "../_utils/vtkOperations.hpp"

//----------------------------------------------------------------------------
vtkOGSDerivatives::vtkOGSDerivatives() {
	this->field     = NULL;
	this->grad_type = 0;
	this->nProcs    = 0;
	this->procId    = 0;

	this->ComputeDivergence = false;
	this->ComputeCurl       = false;
	this->ComputeQ          = false;

	this->isReqInfo = false;

	#ifdef PARAVIEW_USE_MPI
		this->Controller = NULL;
		this->SetController(vtkMultiProcessController::GetGlobalController());
	#endif
}

//----------------------------------------------------------------------------
vtkOGSDerivatives::~vtkOGSDerivatives() {
	this->Setfield(NULL);
	
	#ifdef PARAVIEW_USE_MPI
		this->SetController(NULL);	
	#endif
}

//----------------------------------------------------------------------------
int vtkOGSDerivatives::RequestData( vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector, vtkInformationVector *outputVector) {
	// Get the info objects
	vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
	vtkInformation *outInfo = outputVector->GetInformationObject(0);

	// Stop all threads except from the master to execute
	#ifdef PARAVIEW_USE_MPI
	if (this->procId > 0) return 1;
	#endif

	// Get the input and output
	vtkRectilinearGrid *input = vtkRectilinearGrid::SafeDownCast(
		inInfo->Get(vtkDataObject::DATA_OBJECT()));
	vtkRectilinearGrid *output = vtkRectilinearGrid::SafeDownCast(
		outInfo->Get(vtkDataObject::DATA_OBJECT()));

	output->ShallowCopy(input);

	this->UpdateProgress(0.);

	// Try to load the array, if we fail it is a point array and we should
	// inform the user that point arrays are not treated in this function
	VTKARRAY *vtkArray;
	vtkArray = VTKARRAY::SafeDownCast(input->GetCellData()->GetArray(this->field));
	if (vtkArray == NULL) {
		vtkErrorMacro("Input array "<<this->field<<"is not a cell array."<<
			"This function can only deal with cell arrays!");
		return 0;
	}

	// Here we have a valid vtkArray. Convert it to a field.
	field::Field<FLDARRAY> array;
	array = VTK::createFieldfromVTK<VTKARRAY,FLDARRAY>(vtkArray);
	std::string arrName = std::string(vtkArray->GetName());
	
	// Also, recover XYZ number of nodes
	int nx = input->GetXCoordinates()->GetNumberOfTuples();
	int ny = input->GetYCoordinates()->GetNumberOfTuples();
	int nz = input->GetZCoordinates()->GetNumberOfTuples();

	// Also, recover the weights
	VTKARRAY *vtke1 = NULL, *vtke2 = NULL, *vtke3 = NULL;

	vtke1 = VTKARRAY::SafeDownCast(input->GetCellData()->GetArray("e1"));
	vtke2 = VTKARRAY::SafeDownCast(input->GetCellData()->GetArray("e2"));
	vtke3 = VTKARRAY::SafeDownCast(input->GetCellData()->GetArray("e3"));

	if (this->grad_type > 1 && (vtke1 == NULL || vtke2 == NULL || vtke3 == NULL)) {
		vtkErrorMacro("Mesh weights (e1, e2 and e3) need to be loaded to proceed!");
		return 0;
	}

	// Convert to field arrays
	field::Field<FLDARRAY> e1, e2, e3;
	if (vtke1) e1 = VTK::createFieldfromVTK<VTKARRAY,FLDARRAY>(vtke1);
	if (vtke2) e2 = VTK::createFieldfromVTK<VTKARRAY,FLDARRAY>(vtke2);
	if (vtke3) e3 = VTK::createFieldfromVTK<VTKARRAY,FLDARRAY>(vtke3);

	// Cell centers, update only when request info
	if (this->isReqInfo || this->xyz.isempty()) {
		this->isReqInfo = false;
		
		// Recover Metadata array (depth factor)
		vtkStringArray *vtkmetadata = vtkStringArray::SafeDownCast(
			input->GetFieldData()->GetAbstractArray("Metadata"));
		double dfact = (vtkmetadata != NULL) ? std::stod( vtkmetadata->GetValue(2) ) : 1000.;
		
		if (vtkmetadata == NULL) 
			vtkWarningMacro("Field array Metadata not found! Depth factor set to 1000. automatically.");

		this->xyz = VTK::getVTKCellCenters(input,dfact);
	}

	// Gradient, Divergence, Curl and Q
	field::Field<FLDARRAY> grad, div, curl, Q;

	// Divergence computation
	if (this->ComputeDivergence) {
		if (array.get_m() == 3)
			div.set_dim(array.get_n(),1);
		else
			vtkWarningMacro("Cannot compute the divergence of a scalar array!");
	}
	
	// Curl computation
	if (this->ComputeCurl) {
		if (array.get_m() == 3)
			curl.set_dim(array.get_n(),3);
		else
			vtkWarningMacro("Cannot compute the curl of a scalar array!");
	}

	// Q-criterion computation
	if (this->ComputeQ) {
		if (array.get_m() == 3)
			Q.set_dim(array.get_n(),1);
		else
			vtkWarningMacro("Cannot compute the Q-criterion of a scalar array!");
	}

	// Selection of the gradient method
	switch (this->grad_type) {
		case 0: // Second order, face centered gradient
				// This gradient is unsafe as it relies on the mesh projection
			grad = field::gradXYZ2(nx-1,ny-1,nz-1,this->xyz,array,div,curl,Q);
			break;
		case 1: // Fourth order, face centered gradient
				// This gradient is unsafe as it relies on the mesh projection
			grad = field::gradXYZ4(nx-1,ny-1,nz-1,this->xyz,array,div,curl,Q);
			break;
		case 2:	// OGSTM-BFM approach according to the NEMO handbook
				// This gradient is safe as it relies on the code implementation
			grad = field::gradOGS1(nx-1,ny-1,nz-1,array,e1,e2,e3,div,curl,Q);
			break;
		case 3:	// 2nd order OGSTM-BFM approach
				// This gradient is experimental
			grad = field::gradOGS2(nx-1,ny-1,nz-1,array,e1,e2,e3,div,curl,Q);
			break;
		case 4:	// 4th order OGSTM-BFM approach
				// This gradient is experimental
			grad = field::gradOGS4(nx-1,ny-1,nz-1,array,e1,e2,e3,div,curl,Q);
			break;
		default:
			vtkErrorMacro("Oops! Trouble selecting the gradient type! This should never have happened");
			return 0;
	}

	// Add arrays to output
	std::string gradName = std::string("grad(") + arrName + std::string(")");
	vtkArray = VTK::createVTKfromField<VTKARRAY,FLDARRAY>(gradName.c_str(),grad);
	output->GetCellData()->AddArray(vtkArray);
	vtkArray->Delete();

	// Add Divergence, Curl and Q if computed
	if (!div.isempty()) {
		std::string divName = std::string("div(") + arrName + std::string(")");
		vtkArray = VTK::createVTKfromField<VTKARRAY,FLDARRAY>(divName.c_str(),div);
		output->GetCellData()->AddArray(vtkArray);
		vtkArray->Delete();		
	}

	if (!curl.isempty()) {
		std::string curlName = std::string("curl(") + arrName + std::string(")");
		vtkArray = VTK::createVTKfromField<VTKARRAY,FLDARRAY>(curlName.c_str(),curl);
		output->GetCellData()->AddArray(vtkArray);
		vtkArray->Delete();		
	}

	if (!Q.isempty()) {
		std::string QName = std::string("Q(") + arrName + std::string(")");
		vtkArray = VTK::createVTKfromField<VTKARRAY,FLDARRAY>(QName.c_str(),Q);
		output->GetCellData()->AddArray(vtkArray);
		vtkArray->Delete();		
	}

	// Copy the input grid
	this->UpdateProgress(1.);
	return 1;
}

//----------------------------------------------------------------------------
int vtkOGSDerivatives::RequestInformation(vtkInformation* vtkNotUsed(request),
  vtkInformationVector** vtkNotUsed(inputVector), vtkInformationVector* outputVector) {

	/* SET UP THE PARALLEL CONTROLLER

		The MPI threads come initialized by the ParaView server. Here
		we set up the environment for this filter.

	*/
	#ifdef PARAVIEW_USE_MPI
	if (this->Controller->GetNumberOfProcesses() > 1) {
		this->nProcs = this->Controller->GetNumberOfProcesses();
		this->procId = this->Controller->GetLocalProcessId();
	}

	// Stop all threads except from the master to execute
	if (this->procId > 0) return 1;
	#endif

	this->isReqInfo = true;
	return 1;
}
