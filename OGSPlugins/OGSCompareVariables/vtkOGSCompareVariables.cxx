/*=========================================================================

  Program:   OGSCompareVariables
  Module:    vtkOGSCompareVariables.cxx

  Copyright (c) 2019 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkOGSCompareVariables.h"

#include "vtkDataSet.h"
#include "vtkCellData.h"
#include "vtkPointData.h"
#include "vtkDataArray.h"
#include "vtkFloatArray.h"
#include "vtkDoubleArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"

#include <string>
#include <sstream>
#include <vector>

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
vtkCxxSetObjectMacro(vtkOGSCompareVariables, Controller, vtkMultiProcessController);
#endif

vtkStandardNewMacro(vtkOGSCompareVariables);

//----------------------------------------------------------------------------

/*
	Macro to set the array precision 
*/
#define FLDARRAY double
#define VTKARRAY vtkDoubleArray

#include "../_utils/field.h"
#include "../_utils/vtkFields.hpp"

//----------------------------------------------------------------------------
std::string& ltrim(std::string& str, const std::string& chars = "\t\n\v\f\r ") { str.erase(0, str.find_first_not_of(chars)); return str; }
std::string& rtrim(std::string& str, const std::string& chars = "\t\n\v\f\r ") { str.erase(str.find_last_not_of(chars) + 1); return str; }
std::string&  trim(std::string& str, const std::string& chars = "\t\n\v\f\r ") { return ltrim(rtrim(str, chars), chars); }

//----------------------------------------------------------------------------
int ParseVariables(char *str, std::vector<std::string> &vars, std::vector<double> &vals) {
	std::istringstream varstream(str);
	std::string line;
	while (std::getline(varstream,line)) {
		try {
			// Parse str by =
			size_t pos = line.find("=");
			std::string t1 = line.substr(0, pos);
			std::string t2 = line.substr(pos+1, line.size());
			// Update vectors
			vars.push_back( trim(t1) );
			vals.push_back( std::stod(t2) );
		} catch (...) {
			return 0;
		}
	}
	return 1;
}

//----------------------------------------------------------------------------
void compareVarMax(field::Field<FLDARRAY> &f, std::vector<field::Field<FLDARRAY>> &vars, std::vector<double> &vals, double fillval){
	// Loop the mesh
	#pragma omp parallel
	{
	for (int ii=OMP_THREAD_NUM; ii<f.get_n(); ii+=OMP_NUM_THREADS) {
		// Find the position of the maximum
		double maxval = -1.e20; int posmax = -1;
		for (int jj=0; jj<vars.size(); jj++) {
			posmax = (vars[jj][ii][0] > maxval) ? jj : posmax;
			maxval = (vars[jj][ii][0] > maxval) ? vars[jj][ii][0] : maxval;
		}
		// Fill the array
		f[ii][0] = (posmax >= 0) ? vals[posmax] : fillval;
	}
	}	
}

//----------------------------------------------------------------------------
void compareVarMin(field::Field<FLDARRAY> &f, std::vector<field::Field<FLDARRAY>> &vars, std::vector<double> &vals, double fillval){
	// Loop the mesh
	#pragma omp parallel
	{
	for (int ii=OMP_THREAD_NUM; ii<f.get_n(); ii+=OMP_NUM_THREADS) {
		// Find the position of the maximum
		double minval = 1.e20; int posmin = -1;
		for (int jj=0; jj<vars.size(); jj++) {
			posmin = (vars[jj][ii][0] < minval) ? jj : posmin;
			minval = (vars[jj][ii][0] < minval) ? vars[jj][ii][0] : minval;
		}
		// Fill the array
		f[ii][0] = (posmin >= 0) ? vals[posmin] : fillval;
	}
	}	
}

//----------------------------------------------------------------------------
vtkOGSCompareVariables::vtkOGSCompareVariables() {
	this->variables = NULL;
	this->mode      = 0;
	this->nProcs    = 0;
	this->procId    = 0;
	this->defval    = 0.;

	#ifdef PARAVIEW_USE_MPI
		this->Controller = NULL;
		this->SetController(vtkMultiProcessController::GetGlobalController());
	#endif
}

//----------------------------------------------------------------------------
vtkOGSCompareVariables::~vtkOGSCompareVariables() {
	this->Setvariables(0);

	#ifdef PARAVIEW_USE_MPI
		this->SetController(NULL);	
	#endif
}

//----------------------------------------------------------------------------
int vtkOGSCompareVariables::RequestInformation(vtkInformation* vtkNotUsed(request),
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

	return 1;
}

//----------------------------------------------------------------------------
int vtkOGSCompareVariables::RequestData(vtkInformation *vtkNotUsed(request), 
	vtkInformationVector **inputVector, vtkInformationVector *outputVector) {

	// get the info objects
	vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
	vtkInformation *outInfo = outputVector->GetInformationObject(0);

	// Stop all threads except from the master to execute
	#ifdef PARAVIEW_USE_MPI
	if (this->procId > 0) return 1;
	#endif

	// get the input and output
	vtkDataSet *input = vtkDataSet::SafeDownCast(
		inInfo->Get(vtkDataObject::DATA_OBJECT()));
	vtkDataSet *output = vtkDataSet::SafeDownCast(
		outInfo->Get(vtkDataObject::DATA_OBJECT()));

	output->ShallowCopy(input);
	this->UpdateProgress(0.);

	// Parse user inputted variables
	std::vector<std::string> vars;
	std::vector<double> vals;
	if (!ParseVariables(this->variables,vars,vals)) {
		vtkErrorMacro("Error parsing variables dialog. Please review synthax.");
		return 0;
	}
	this->UpdateProgress(0.25);

	// Try to load the variables the user has inputted. If the variable is not
	// recognized stop and throw an error.
	std::vector<field::Field<FLDARRAY>> varfield;
	VTKARRAY *vtkArray = NULL; bool celldata = true;
	for(std::string v : vars) {
		vtkArray = VTKARRAY::SafeDownCast(input->GetCellData()->GetArray( v.c_str() ));
		if (!vtkArray) {
			vtkArray = VTKARRAY::SafeDownCast(input->GetPointData()->GetArray( v.c_str() ));
			celldata = false;
		}
		// We couldn't find the array
		if (!vtkArray) { 
			vtkErrorMacro("Could not find variable <" << v <<">");
			return 0;
		}
		// We do have the array here, so update the vector
		varfield.push_back( VTK::createFieldfromVTK<VTKARRAY,FLDARRAY>(vtkArray) );
	}
	this->UpdateProgress(0.5);

	// Define an array to contain the comparison
	field::Field<FLDARRAY> comp(varfield[0].get_n(),1,this->defval);

	// Get the maximum or minimum of an array according to the mode of operation
	switch (this->mode) {
		case 0:
			compareVarMax(comp,varfield,vals,this->defval);
			break;
		case 1:
			compareVarMin(comp,varfield,vals,this->defval);
			break;
	}
	this->UpdateProgress(0.75);

	// Add array to output
	vtkArray = VTK::createVTKfromField<VTKARRAY,FLDARRAY>("comp",comp);
	// Add the new variable to the input
	if (celldata)
		output->GetCellData()->AddArray(vtkArray);
	else
		output->GetPointData()->AddArray(vtkArray);

	// Update progress and leave
	this->UpdateProgress(1.);
	return 1;
}
