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
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"

#include <string>
#include <sstream>
#include <algorithm>
#include <vector>

#include "pugixml.hpp"
namespace xml = pugi;

vtkStandardNewMacro(vtkOGSVariableAggregator);

class vtkOGSVariableAggregator::vtkVectorOfArrays :
  public std::vector<vtkFloatArray*>
{
};

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
}

//----------------------------------------------------------------------------
vtkOGSVariableAggregator::~vtkOGSVariableAggregator() {
	this->VarDataArraySelection->Delete();

	this->SetFileName(0);
	this->SetXMLText(0);
}

//----------------------------------------------------------------------------
int vtkOGSVariableAggregator::RequestData(vtkInformation *vtkNotUsed(request), 
	vtkInformationVector **inputVector, vtkInformationVector *outputVector) {

	// get the info objects
	vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
	vtkInformation *outInfo = outputVector->GetInformationObject(0);

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
	for(int varId=0; varId<this->GetNumberOfVarArrays(); varId++) {
		// Check if the variable has been enabled
		if ( !this->GetVarArrayStatus(this->GetVarArrayName(varId)) ) continue;
		// Create a new vtkFloatArray to store the new variable
		vtkFloatArray *vtkArray = vtkFloatArray::New();
		vtkArray->SetName(this->GetVarArrayName(varId));
		vtkArray->SetNumberOfComponents(1); // Scalar field
		// Obtain the variables to aggregate
		std::string vararray = this->AggrVar[this->GetVarArrayName(varId)];
		int celldata = 1;
		int current,previous=0;
		vtkVectorOfArrays *AgrVarArray;
		AgrVarArray = new vtkVectorOfArrays;
		do {
			// Find an occurence of the delimiter and split the string
			current = vararray.find(';',previous);
			std::string varname = vararray.substr(previous, current-previous);
			if (varname == "") continue;
			// Now varname contains the name of the variable to load
			// Try to load the array as cell data
			vtkDataArray *array;
			array = input->GetCellData()->GetArray(varname.c_str());
			if (!array) {// Then array might be point data
				array    = input->GetPointData()->GetArray(varname.c_str());
				celldata = 0;
			}
			if (!array) { // We couldn't find the array
				vtkErrorMacro("Could not find variable <" << varname <<">");
				return 0;
			}
			// Array should exist at this point, store it
			AgrVarArray->push_back( vtkFloatArray::SafeDownCast(array) );
			previous = current + 1;
		}while(current != std::string::npos);
		// Here we have the array to aggregate and the arrays to aggregate from
		// Get the size of the mesh
		size_t meshsize = AgrVarArray->at(0)->GetNumberOfTuples();
		vtkArray->SetNumberOfTuples(meshsize);
		// Loop the mesh and create the new variable
		for (int ii=0; ii<meshsize; ii++) {
			double aux = 0.;
			vtkVectorOfArrays::iterator iter;
			for (iter=AgrVarArray->begin(); iter != AgrVarArray->end(); iter++)
				aux += (*iter)->GetTuple1(ii);
			vtkArray->SetTuple1(ii,aux);
		}
		// Add the new variable to the input
		if (celldata)
			output->GetCellData()->AddArray(vtkArray);
		else
			output->GetPointData()->AddArray(vtkArray);
		// Update progress
		this->UpdateProgress(1./(double)(this->GetNumberOfVarArrays())*varId);
		// Delete
		vtkArray->Delete();
		delete AgrVarArray;
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
int vtkOGSVariableAggregator::RequestInformation(vtkInformation* vtkNotUsed(request),
  vtkInformationVector** vtkNotUsed(inputVector), vtkInformationVector* outputVector) {
	// Parse XML and TextBox
	this->ParseXML();
	this->SetAggrVarsText();

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
	if ( doc.load_string(this->XMLText) ) {
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