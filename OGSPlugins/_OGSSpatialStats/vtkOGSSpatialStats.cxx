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


namespace VTK
{
// Include the VTK Operations
#include "../_utils/vtkOperations.h"
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

	// Compute the cell centers and number of different z points
//	int nz = 0;
//	std::vector<double> zcoords;
//	vtkFloatArray *vtkCellCenters;
//
//	vtkCellCenters = ComputeCellCenters(input,this->depth_factor,nz,zcoords);

	// Loop the number of variables
	// For each variable, we will see if a prost processing exists and
	// then we will loop the mesh and create the statistics.
//	int narrays = input->GetCellData()->GetNumberOfArrays();
//	int nstat   = this->GetNumberOfStatArrays();
//	for (int varId = 0; varId < narrays; varId++) {
//		// Recover the array and the array name
//		vtkFloatArray *vtkVarArray = vtkFloatArray::SafeDownCast(
//			input->GetCellData()->GetArray(varId));
//		char *array_name = vtkVarArray->GetName();
//
//		// We need to decide which statistics to compute. Create a map and initialize
//		// a map that will link the type of statistic and the associated vtkFloatArray
//		std::map<std::string, vtkFloatArray*> mapStatArray;
//
//		for (int statId = 0; statId < nstat; statId++) {
//			const char *statName = this->GetStatArrayName(statId);
//			// Skip those arrays who have not been enabled
//			if (!this->GetStatArrayStatus(statName))
//				continue;
//			// Statistical variable name
//			char statVarName[256];
//			sprintf(statVarName,"%s, %s",array_name,statName);
//			// Define a new vtkFloatArray
//			vtkFloatArray *vtkStatVar = vtkFloatArray::New();
//			vtkStatVar->SetName(statVarName);
//			vtkStatVar->SetNumberOfComponents( vtkVarArray->GetNumberOfComponents() );
//			vtkStatVar->SetNumberOfTuples(ncells);
//			vtkStatVar->Fill(0.);
//			// Store the array in the map
//			mapStatArray.insert(std::make_pair(std::string(statName),vtkStatVar));
//		}
//
//		// Define and allocate parameters
//		int *nz_dep;
//		double *meanval, *maxval, *minval;
//		meanval = (double*)malloc(nz*sizeof(double));
//		nz_dep  = (int*)malloc(nz*sizeof(int));
//		maxval  = (double*)malloc(nz*sizeof(double));
//		minval  = (double*)malloc(nz*sizeof(double));
//		
//		for (int kk = 0; kk < nz; kk++) {
//			meanval[kk] = 0.;
//			maxval[kk]  = -1.e20;
//			minval[kk]  = 1.e20;
//		}
//
//		// First loop in the mesh
//		// Compute meanval, maxval and minval per depth level
//		for (int cellId = 0; cellId < ncells; cellId++) {
//			// Recover the cell center
//			double xyz[3];
//			vtkCellCenters->GetTuple(cellId,xyz);
//			// Recover the variable value
//			double value = vtkVarArray->GetTuple1(cellId);
//			// Obtain the id for the current depth
//			int zId = -1;
//			for (int kk = 0; kk < nz && zId < 0; kk++)
//				if ( fabs(xyz[2] - zcoords.at(kk)) < 1.e-4 ) zId = kk;
//			if (zId < 0)
//				vtkErrorMacro("Couldn't find depth level");
//			// Set the variable according to the depth level
//			// Mean
//			meanval[zId] += value;
//			nz_dep[zId]++;
//			// Maximum and minimum
//			if (value > maxval[zId]) maxval[zId] = value;
//			if (value < minval[zId]) minval[zId] = value;
//		}
//		for (int kk = 0; kk < nz; kk++) 
//			meanval[kk] /= (double)(nz_dep[kk]);
//
//
//		// Final mesh loop
//		// Set values to array
//		for (int cellId = 0; cellId < ncells; cellId++) {
//			// Recover the cell center
//			double xyz[3];
//			vtkCellCenters->GetTuple(cellId,xyz);
//			// Recover the variable value
//			double value = vtkVarArray->GetTuple1(cellId);
//			// Obtain the id for the current depth
//			int zId = -1;
//			for (int kk = 0; kk < nz && zId < 0; kk++)
//				if ( fabs(xyz[2] - zcoords.at(kk)) < 1.e-4 ) zId = kk;
//			if (zId < 0)
//				vtkErrorMacro("Couldn't find depth level");
//			// Set values for statistics
//			std::map<std::string,vtkFloatArray*>::iterator iter;
//			for (iter = mapStatArray.begin(); iter != mapStatArray.end(); iter++) {
//				// Which statistic are we computing?
//				int statId = this->GetStatArrayIndex(iter->first.c_str());
//				switch (statId) {
//					case 0: // Mean
//						iter->second->SetTuple1(cellId,meanval[zId]);
//						break;
//					case 1: // Std dev
//						break;
//					case 2: // Min
//						iter->second->SetTuple1(cellId,minval[zId]);
//						break;
//					case 3: // p05
//						break;
//					case 4: // p25
//						break;
//					case 5: // p50
//						break;
//					case 6: // p75
//						break;
//					case 7: // p95
//						break;
//					case 8: // Max
//						iter->second->SetTuple1(cellId,maxval[zId]);
//						break;
//				}
//			}
//		}
//
//
//
//		free(meanval); free(nz_dep);
//		free(maxval); free(minval);

//		double meanval = 0, maxval = -1.e20, minval = 1.e20;
//		for (int cellId = 0; cellId < ncells; cellId++) {
//			double value = vtkVarArray->GetTuple1(cellId);
//			// Mean
//			meanval += value;
//			// Max and min
//			maxval = value > maxval ? value : maxval;
//			minval = value < minval ? value : minval;
//		}
//		meanval /= (double)(ncells);

		// Second mesh loop
		// Compute stddev and set values
//		for (int cellId = 0; cellId < ncells; cellId++) {
//			// stddev
//			// Set values for statistics
//			std::map<std::string,vtkFloatArray*>::iterator iter;
//			for (iter = mapStatArray.begin(); iter != mapStatArray.end(); iter++) {
//				// Which statistic are we computing?
//				int statId = this->GetStatArrayIndex(iter->first.c_str());
//				switch (statId) {
//					case 0: // Mean
//						iter->second->SetTuple1(cellId,meanval);
//						break;
//					case 1: // Std dev
//						break;
//					case 2: // Min
//						iter->second->SetTuple1(cellId,minval);
//						break;
//					case 3: // p05
//						break;
//					case 4: // p25
//						break;
//					case 5: // p50
//						break;
//					case 6: // p75
//						break;
//					case 7: // p95
//						break;
//					case 8: // Max
//						iter->second->SetTuple1(cellId,maxval);
//						break;
//				}
//			}
//		}

//		// Now that we computed the arrays, we can set them in the output
//		// and deallocate memory
//		std::map<std::string,vtkFloatArray*>::iterator iter;
//		for (iter = mapStatArray.begin(); iter != mapStatArray.end(); iter++) {
//			output->GetCellData()->AddArray(iter->second);
//			iter->second->Delete();
//		}
//
//		this->UpdateProgress(0.+1./(double)(narrays)*(double)(varId));
//	}
	// Deallocate cell centers
	vtkCellCenters->Delete();
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
