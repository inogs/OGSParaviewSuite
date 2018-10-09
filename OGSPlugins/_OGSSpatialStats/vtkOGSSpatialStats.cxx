/*=========================================================================

  Program:   OGSSpatialStats
  Module:    vtkOGSSpatialStats.cxx

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkFloatArray.h"
#include "vtkStringArray.h"
#include "vtkCell.h"
#include "vtkPoints.h"
#include "vtkCellData.h"
#include "vtkFieldData.h"
#include "vtkDataSet.h"
#include "vtkRectilinearGrid.h"
#include "vtkDataArraySelection.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#include "vtkOGSSpatialStats.h"

#include "vtkObjectFactory.h"

#include <string>
#include <algorithm>
#include <map>
#include <vector>
#include <numeric>


namespace VTK
{
// Include the VTK Operations
#include "../_utils/vtkOperations.cpp"
}


vtkStandardNewMacro(vtkOGSSpatialStats);

//----------------------------------------------------------------------------
vtkOGSSpatialStats::vtkOGSSpatialStats(){
	this->StatDataArraySelection = vtkDataArraySelection::New();
	this->StatDataArraySelection->AddArray("Mean");
	this->StatDataArraySelection->AddArray("std");
	this->StatDataArraySelection->AddArray("min");
	this->StatDataArraySelection->AddArray("p05");
	this->StatDataArraySelection->AddArray("p25");
	this->StatDataArraySelection->AddArray("p50");
	this->StatDataArraySelection->AddArray("p75");
	this->StatDataArraySelection->AddArray("p95");
	this->StatDataArraySelection->AddArray("max");

	this->depth_factor = 1000.;
}

//----------------------------------------------------------------------------
vtkOGSSpatialStats::~vtkOGSSpatialStats() {
	this->StatDataArraySelection->Delete();
}

//----------------------------------------------------------------------------
int vtkOGSSpatialStats::RequestData(vtkInformation *vtkNotUsed(request),
	vtkInformationVector **inputVector,vtkInformationVector *outputVector) {
	// Get the info objects
	vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
	vtkInformation *outInfo = outputVector->GetInformationObject(0);

	// Get the input and output
	vtkDataSet *input = vtkDataSet::SafeDownCast(
		inInfo->Get(vtkDataObject::DATA_OBJECT()));
	vtkDataSet *output = vtkDataSet::SafeDownCast(
		outInfo->Get(vtkDataObject::DATA_OBJECT()));

	// We just want to copy the mesh, not the variables
	output->CopyStructure(input); 

	this->UpdateProgress(0.);

	// At this point we either have a rectilinear grid or an
	// unstructured grid with either cell or point arrays

	// Start by dealing with any cell array
	this->CellStats(input,output,0.0);
	this->UpdateProgress(1.);

	// Start by dealing with any cell array
//	this->PointStats(input,output,0.5);
	this->UpdateProgress(1.);

	return 1;
}

//----------------------------------------------------------------------------
void vtkOGSSpatialStats::CellStats(vtkDataSet *input, vtkDataSet *output, double updst) {
	/*
		Computation of spatial statistics on cell data.

		Grabs the input and returns an output with the arrays filled.
	*/
	int ncells = input->GetNumberOfCells();

	// Compute the cell centers
	vtkFloatArray *vtkCellCenters = VTK::getCellCoordinates("Cell centers",input);
	this->UpdateProgress(updst+0.05);

	// Count the number of unique z coordinates
	int *cellId2zId; cellId2zId = (int*)malloc(ncells*sizeof(int));
	std::vector<double> zcoords;
// TODO: epsilon should not be hardcoded
// TODO: implement range of depths
	int nz = VTK::countUniqueZ(vtkCellCenters,ncells,1000.*1e-3,cellId2zId,zcoords);
	this->UpdateProgress(updst+0.1);

	// The previous operations could be grouped into one mesh loop
	// For the sake of simplicity and reusability, they have been separated
	// Moreover, this proved to be a faster approach.

	// Recover e1t and e2t
	vtkFloatArray *vtke1t = vtkFloatArray::SafeDownCast(
		input->GetCellData()->GetArray("e1t"));
	output->GetCellData()->AddArray(vtke1t);
	vtkFloatArray *vtke2t = vtkFloatArray::SafeDownCast(
		input->GetCellData()->GetArray("e2t"));
	output->GetCellData()->AddArray(vtke2t);

	// We can now loop the number of active variables
	int narrays = input->GetCellData()->GetNumberOfArrays();
	int nstat   = this->GetNumberOfStatArrays();

	for (int varId = 0; varId < narrays; varId++) {
		// Recover the array and the array name
		vtkFloatArray *vtkVarArray = vtkFloatArray::SafeDownCast(
			input->GetCellData()->GetArray(varId));
		char *array_name = vtkVarArray->GetName();

		// We should not average the coast or basins mask nor e1t or e2t
		// Names have been harcoded here as there is no way to ensure that
		// these arrays will exist or not.
		if (std::string(array_name) == "coast mask")  continue;
		if (std::string(array_name) == "basins mask") continue;
		if (std::string(array_name) == "e1t")         continue;
		if (std::string(array_name) == "e2t")         continue;

		// Also anything not being an scalar array should not be computed
		if (vtkVarArray->GetNumberOfComponents() > 1) {
			vtkWarningMacro("Variable "<<array_name<<
				" is not an scalar array. Statistics will be skipped for this variable.");
			continue;
		}

		// We need to decide which statistics to compute. Create a map and initialize
		// a map that will link the type of statistic and the associated vtkFloatArray
		std::map<std::string, vtkFloatArray*> mapStatArray;

		for (int statId = 0; statId < nstat; statId++) {
			const char *statName = this->GetStatArrayName(statId);
			// Skip those arrays who have not been enabled
			if (!this->GetStatArrayStatus(statName))
				continue;
			// Statistical variable name
			char statVarName[256];
			sprintf(statVarName,"%s, %s",array_name,statName);
			// Define a new vtkFloatArray
			vtkFloatArray *vtkStatVar = vtkFloatArray::New();
			vtkStatVar->SetName(statVarName);
			vtkStatVar->SetNumberOfComponents( vtkVarArray->GetNumberOfComponents() );
			vtkStatVar->SetNumberOfTuples(ncells);
			vtkStatVar->Fill(0.);
			// Store the array in the map
			mapStatArray.insert(std::make_pair(std::string(statName),vtkStatVar));
		}

		// Create a map that where for each depth level it will store the values
		// of the variables and the weight as a 2D layer
		std::map<int,std::vector<double>> vPerLayer;
		std::map<int,std::vector<double>> wPerLayer;
		std::map<int,std::vector<int>>    cPerLayer; // Cell per layer

		// First loop in the mesh
		// Set the vPerLayer and wPerLayer
		for (int cellId = 0; cellId < ncells; cellId++) {
			// Recover the variable value
			double value = vtkVarArray->GetTuple1(cellId);
			double e1t   = vtke1t->GetTuple1(cellId);
			double e2t   = vtke2t->GetTuple1(cellId);
			// Obtain the id for the current depth
			int zId = cellId2zId[cellId];
			if (zId < 0) vtkErrorMacro("Error computing <cellId2zId>");
			// Store the variable per each layer
			vPerLayer[zId].push_back(value);
			wPerLayer[zId].push_back(e1t*e2t);
			cPerLayer[zId].push_back(cellId);
		}

		// Loop the depth layers
		for (int kk = 0; kk < nz; kk++) {
			// Number of elements per layer
			int nlayer = vPerLayer[kk].size();
			// Occurrences vector
			std::vector< std::pair<double,int> > orderVal;

			// First loop on the layer, compute the weights, mean and min/max
			double sum_weight = 0., meanval = 0., maxval = -1.e20, minval = 1.e20;
			for (int ii = 0; ii < nlayer; ii++) {
				// Minimum and Maximum
				minval = (vPerLayer[kk].at(ii) < minval) ? vPerLayer[kk].at(ii) : minval;
				maxval = (vPerLayer[kk].at(ii) > maxval) ? vPerLayer[kk].at(ii) : maxval;
				// Weights and Mean
				sum_weight += wPerLayer[kk].at(ii);
				meanval += vPerLayer[kk].at(ii)*wPerLayer[kk].at(ii);
				// List to order values
				orderVal.push_back( std::make_pair(vPerLayer[kk].at(ii),ii) );
			}
			// Mean
			meanval /= sum_weight;

			// Order the values
			std::sort(orderVal.begin(),orderVal.end());

			// Second loop on the layer, std and percentile weight
			double stdval = 0., weights = 0.;
			std::vector<double> percw;
			for (int ii = 0; ii < nlayer; ii++) {
				// Std
				double aux = vPerLayer[kk].at(ii) - meanval;
				stdval += aux*aux*wPerLayer[kk].at(ii);
				// Percentile weight
				int ind = orderVal[ii].second;
				weights += wPerLayer[kk].at(ind);
				aux = (weights-0.5*wPerLayer[kk].at(ind))/sum_weight; // Reused variable
				percw.push_back( aux );
			}
			// Std
			stdval = sqrt(stdval/sum_weight);

			// Compute the percentiles
			double perc[]    = {.05,.25,.50,.75,.95};
			double percval[] = { 0., 0., 0., 0., 0.};

			for (int pp = 0; pp < 5; pp++) {
				// Find the value that is equal to perc or immediately after.
				std::vector<double>::iterator lbound;
				lbound = std::lower_bound(percw.begin(),percw.end(),perc[pp]);
				lbound--; // We need to decrement this value;
				int s = lbound - percw.begin(); // This is our position on the ordered value array

				// Set the value for the weight
				if (s == 0)        {percval[pp] = orderVal[s].first; continue;} // == sd[0]
				if (s == nlayer-1) {percval[pp] = orderVal[s].first; continue;} // == sd[n-1]

				double f1 = (percw[s] - perc[pp])   / (percw[s] - percw[s-1]);
				double f2 = (perc[pp] - percw[s-1]) / (percw[s] - percw[s-1]);

				percval[pp] = f1*orderVal[s-1].first + f2*orderVal[s].first;
			}

			// Third loop on the layer, set the mesh
			for (int ii = 0; ii < nlayer; ii++) {
				int cellId = cPerLayer[kk].at(ii);
				// Set values for statistics
				std::map<std::string,vtkFloatArray*>::iterator iter;
				for (iter = mapStatArray.begin(); iter != mapStatArray.end(); iter++) {
					// Which statistic are we computing?
					int statId = this->GetStatArrayIndex(iter->first.c_str());
					switch (statId) {
						case 0: // Mean
							iter->second->SetTuple1(cellId,meanval);
							break;
						case 1: // Std dev
							iter->second->SetTuple1(cellId,stdval);
							break;
						case 2: // Min
							iter->second->SetTuple1(cellId,minval);
							break;
						case 3: // p05
							iter->second->SetTuple1(cellId,percval[0]);
							break;
						case 4: // p25
							iter->second->SetTuple1(cellId,percval[1]);
							break;
						case 5: // p50
							iter->second->SetTuple1(cellId,percval[2]);
							break;
						case 6: // p75
							iter->second->SetTuple1(cellId,percval[3]);
							break;
						case 7: // p95
							iter->second->SetTuple1(cellId,percval[4]);
							break;
						case 8: // Max
							iter->second->SetTuple1(cellId,maxval);
							break;
					}
				}
			}
		}

		// Now that we computed the arrays, we can set them in the output
		// and deallocate memory
		std::map<std::string,vtkFloatArray*>::iterator iter;
		for (iter = mapStatArray.begin(); iter != mapStatArray.end(); iter++) {
			output->GetCellData()->AddArray(iter->second);
			iter->second->Delete();
		}

		this->UpdateProgress(updst+0.1+0.4/(double)(narrays)*(double)(varId));
	}
	// Deallocate and delete
	vtkCellCenters->Delete(); free(cellId2zId);
}

//----------------------------------------------------------------------------
void vtkOGSSpatialStats::DisableAllStatArrays()
{
	this->StatDataArraySelection->DisableAllArrays();
}

void vtkOGSSpatialStats::EnableAllStatArrays()
{
	this->StatDataArraySelection->EnableAllArrays();
}

int vtkOGSSpatialStats::GetNumberOfStatArrays()
{
	return this->StatDataArraySelection->GetNumberOfArrays();
}

const char* vtkOGSSpatialStats::GetStatArrayName(int index)
{
	if (index >= (int)this->GetNumberOfStatArrays() || index < 0)
		return NULL;
	else
		return this->StatDataArraySelection->GetArrayName(index);
}

int vtkOGSSpatialStats::GetStatArrayIndex(const char* name)
{
	return this->StatDataArraySelection->GetArrayIndex(name);
}

int vtkOGSSpatialStats::GetStatArrayStatus(const char* name)
{
	return this->StatDataArraySelection->ArrayIsEnabled(name);
}

void vtkOGSSpatialStats::SetStatArrayStatus(const char* name, int status)
{
	if (status)
		this->StatDataArraySelection->EnableArray(name);
	else
		this->StatDataArraySelection->DisableArray(name);

	this->Modified();
}
