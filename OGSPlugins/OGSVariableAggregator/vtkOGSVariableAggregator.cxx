/*=========================================================================

  Program:   OGSVariableAggregator
  Module:    vtkOGSVariableAggregator.cxx

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkOGSVariableAggregator.h"

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

#ifdef PARAVIEW_USE_MPI
#include "vtkMultiProcessController.h"
vtkCxxSetObjectMacro(vtkOGSVariableAggregator, Controller, vtkMultiProcessController);
#endif

vtkStandardNewMacro(vtkOGSVariableAggregator);

//----------------------------------------------------------------------------
#include "macros.h"
#include "field.h"
#include "vtkFields.h"

#include "pugixml.h"
namespace xml = pugi;

//----------------------------------------------------------------------------
void addVar(vtkDataArraySelection *AggrVar, std::map<std::string, std::string> *Var2Aggr) {
	// P_c
	AggrVar->AddArray("P_c");
	Var2Aggr->insert(std::make_pair("P_c","P1c;P2c;P3c;P4c;"));
	// P_l
	AggrVar->AddArray("P_l");
	Var2Aggr->insert(std::make_pair("P_l","P1l;P2l;P3l;P4l;"));
	// T_c
	AggrVar->AddArray("T_c");
	Var2Aggr->insert(std::make_pair("T_c","B1c;P1c;P2c;P3c;P4c;Z3c;Z4c;Z5c;Z6c;"));
	// RTc
	AggrVar->AddArray("RTc");
	Var2Aggr->insert(std::make_pair("RTc","resMEZ1c;resMEZ2c;resMIZ1c;resMIZ2c;resPBAc;resPPY1c;resPPY2c;resPPY3c;resPPY4c;"));
}

//----------------------------------------------------------------------------
vtkOGSVariableAggregator::vtkOGSVariableAggregator() {
	// Define data arrays
	this->VarDataArraySelection = vtkDataArraySelection::New();
	// Preallocate some most used variables
	addVar(this->VarDataArraySelection,&this->AggrVar);

	this->deleteVars = 0;
	this->FileName   = NULL; 
	this->XMLText    = NULL;
	this->nProcs     = 0;
	this->procId     = 0;

	#ifdef PARAVIEW_USE_MPI
		this->Controller = NULL;
		this->SetController(vtkMultiProcessController::GetGlobalController());
	#endif
}

//----------------------------------------------------------------------------
vtkOGSVariableAggregator::~vtkOGSVariableAggregator() {
	this->VarDataArraySelection->Delete();

	this->SetFileName(0);
	this->SetXMLText(0);

	#ifdef PARAVIEW_USE_MPI
		this->SetController(NULL);	
	#endif
}

//----------------------------------------------------------------------------
int vtkOGSVariableAggregator::RequestInformation(vtkInformation* vtkNotUsed(request),
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

	// Parse XML and TextBox
	this->ParseXML();
	this->SetAggrVarsText();

	return 1;
}

//----------------------------------------------------------------------------
int vtkOGSVariableAggregator::RequestData(vtkInformation *vtkNotUsed(request), 
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

	/*
		For each variable in VarDataArraySelection, parse the composite variables
		in ArgDataArraySelection and try to load the variables from the input. Then loop
		the mesh and add the variables at each point. Assume all variables are scalar 
		arrays.
	*/
	// Parallelization strategy MPI
	for(int varId=0; varId<this->GetNumberOfVarArrays(); varId++) {
		// Check if the variable has been enabled
		if ( !this->GetVarArrayStatus(this->GetVarArrayName(varId)) ) continue;

		// Obtain the variables to aggregate
		VTKARRAY *vtkArray = NULL;
		std::string vararray = this->AggrVar[this->GetVarArrayName(varId)];
		int celldata = 1, current = 0, previous = 0;
		std::vector<field::Field<FLDARRAY>> arrayVector;
		do {
			// Find an occurence of the delimiter and split the string
			current = vararray.find(';',previous);
			std::string varname = vararray.substr(previous, current-previous);
			if (varname == "") continue;
			// Now varname contains the name of the variable to load
			// Try to load the array as cell data
			vtkArray = VTKARRAY::SafeDownCast(input->GetCellData()->GetArray(varname.c_str()));
			// Then array might be point data
			if (!vtkArray) { 
				vtkArray = VTKARRAY::SafeDownCast(input->GetPointData()->GetArray(varname.c_str()));
				celldata = 0;
			}
			// We couldn't find the array
			if (!vtkArray) { 
				vtkErrorMacro("Could not find variable <" << varname <<">");
				return 0;
			}
			// Array should exist at this point, store it
			arrayVector.push_back( VTK::createFieldfromVTK<VTKARRAY,FLDARRAY>(vtkArray) );
			// Check we are dealing with scalar arrays
			if (arrayVector.back().get_m() != 1) {
				vtkErrorMacro("Variable <" << varname << "> is not a scalar array. Cannot proceed.");
				return 0;				
			}
			previous = current + 1;
		}while(current != std::string::npos);

		// Here we have the array to aggregate and the arrays to aggregate from
		// then create a new Field to store the new variable
		field::Field<FLDARRAY> arrayNew(arrayVector[0]);
		
		// Loop the vectors and create the new variable (parallelization MPI)
		field::Field<FLDARRAY>::iterator arrIter, auxIter;
		std::vector<field::Field<FLDARRAY>>::iterator vecIter;
		for (vecIter = arrayVector.begin()+1; vecIter != arrayVector.end(); ++vecIter) {
			// Loop the mesh 
			#pragma omp parallel
			{
			for (int ii=OMP_THREAD_NUM; ii<arrayNew.get_n(); ii+=OMP_NUM_THREADS)
				arrayNew[ii][0] += (*vecIter)[ii][0];
			}
		}
		// Convert to vtkArray
		vtkArray = VTK::createVTKfromField<VTKARRAY,FLDARRAY>(this->GetVarArrayName(varId),arrayNew);
		// Add the new variable to the input
		if (celldata) {
			output->GetCellData()->AddArray(vtkArray);
			output->GetCellData()->SetActiveScalars(this->GetVarArrayName(varId));
		} else {
			output->GetPointData()->AddArray(vtkArray);
			output->GetPointData()->SetActiveScalars(this->GetVarArrayName(varId));
		}
		// Update progress
		this->UpdateProgress(1./(double)(this->GetNumberOfVarArrays())*varId);
		// Delete
		vtkArray->Delete();
	}

	/*
		At the user's request, the variables that have been aggregated can
		be deleted. This option is activated by default.
	*/
	if (this->deleteVars) {
		for(int varId=0; varId<this->GetNumberOfVarArrays(); varId++) {
			// Obtain the variables to aggregate
			std::string vararray = this->AggrVar[this->GetVarArrayName(varId)];
			int current,previous=0;
			do {
				// Find an occurence of the delimiter and split the string
				current = vararray.find(';',previous);
				std::string varname = vararray.substr(previous, current-previous);
				// Now varname contains the name of the variable to load
				// Remove the variable
				output->GetCellData()->RemoveArray(varname.c_str());
				output->GetPointData()->RemoveArray(varname.c_str());
				previous = current + 1;
			}while(current != std::string::npos);
		}
	}

	// Update progress and leave
	this->UpdateProgress(1.);
	return 1;
}

//----------------------------------------------------------------------------
void vtkOGSVariableAggregator::ParseXML() {
	// Define an XML doc
	xml::xml_document doc;
	// If reading has been done correctly
	if ( doc.load_file(this->FileName) ) {
		// Define our main node
		xml::xml_node aggregate = doc.child("root").child("vars_for_All_Statistics").child("aggregate");
		// Loop the node
		for(xml::xml_node_iterator iter = aggregate.begin(); iter != aggregate.end(); iter++) {
			// Loop the subnode and create the aggregated array
			std::ostringstream buf;
			for(xml::xml_node_iterator iter2 = iter->begin(); iter2 != iter->end(); iter2++)
				buf << iter2->attribute("name").value() << ";";
			// Define the variable
			const char *varname = iter->attribute("name").value();
			// Check whether the variable exists or not
			if (this->GetVarArrayIndex(varname) < 0) { // Variable doesn't exist
				// Add array
				this->VarDataArraySelection->AddArray(varname);
				this->VarDataArraySelection->EnableArray(varname);
			}
			// Insert or overwrite aggregated variables
			this->AggrVar[varname] = buf.str();
		}
	}
}

//----------------------------------------------------------------------------
void vtkOGSVariableAggregator::SetAggrVarsText() {
	// Define an XML doc
	xml::xml_document doc;
	// If reading has been done correctly
	if ( doc.load(this->XMLText) ) {
		// Define our main node
		xml::xml_node aggregate = doc.child("aggregate");
		// Loop the node
		for(xml::xml_node_iterator iter = aggregate.begin(); iter != aggregate.end(); iter++) {
			// Loop the subnode and create the aggregated array
			std::ostringstream buf;
			for(xml::xml_node_iterator iter2 = iter->begin(); iter2 != iter->end(); iter2++)
				buf << iter2->attribute("name").value() << ";";
			// Define the variable
			const char *varname = iter->attribute("name").value();
			// Check whether the variable exists or not
			if (this->GetVarArrayIndex(varname) < 0) { // Variable doesn't exist
				// Add array
				this->VarDataArraySelection->AddArray(varname);
				this->VarDataArraySelection->EnableArray(varname);
			}
			// Insert or overwrite aggregated variables
			this->AggrVar[varname] = buf.str();
		}
	}
}

//----------------------------------------------------------------------------
void vtkOGSVariableAggregator::DisableAllVarArrays()
{
	this->VarDataArraySelection->DisableAllArrays();
}

void vtkOGSVariableAggregator::EnableAllVarArrays()
{
	this->VarDataArraySelection->EnableAllArrays();
}

int vtkOGSVariableAggregator::GetNumberOfVarArrays()
{
	return this->VarDataArraySelection->GetNumberOfArrays();
}

const char* vtkOGSVariableAggregator::GetVarArrayName(int index)
{
	if (index >= (int)this->GetNumberOfVarArrays() || index < 0)
		return NULL;
	else
		return this->VarDataArraySelection->GetArrayName(index);
}

int vtkOGSVariableAggregator::GetVarArrayIndex(const char* name)
{
	return this->VarDataArraySelection->GetArrayIndex(name);
}

int vtkOGSVariableAggregator::GetVarArrayStatus(const char* name)
{
	return this->VarDataArraySelection->ArrayIsEnabled(name);
}

void vtkOGSVariableAggregator::SetVarArrayStatus(const char* name, int status)
{
	if (status)
		this->VarDataArraySelection->EnableArray(name);
	else
		this->VarDataArraySelection->DisableArray(name);

	this->Modified();
}
