/*=========================================================================

  Program:   OGSSpatialStats
  Module:    vtkOGSSpatialStats.cxx

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkIntArray.h"
#include "vtkFloatArray.h"
#include "vtkDoubleArray.h"
#include "vtkStringArray.h"
#include "vtkCell.h"
#include "vtkPoints.h"
#include "vtkCellData.h"
#include "vtkPointData.h"
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
#include <vector>
#include <map>
#include <algorithm>

#ifdef PARAVIEW_USE_MPI
#include "vtkMultiProcessController.h"
vtkCxxSetObjectMacro(vtkOGSSpatialStats, Controller, vtkMultiProcessController);
#endif

vtkStandardNewMacro(vtkOGSSpatialStats);

//----------------------------------------------------------------------------
#include "macros.h"
#include "fieldOperations.h"
#include "vtkFields.h"
#include "vtkOperations.h"

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

	this->iscelld   = true;
	this->isReqInfo = false;
	this->epsi      = 1.e-3;
	this->ndepths   = 0;
	this->useVolume = 0;
	this->nProcs    = 0;
	this->procId    = 0;

	#ifdef PARAVIEW_USE_MPI
		this->Controller = NULL;
		this->SetController(vtkMultiProcessController::GetGlobalController());
	#endif
}

//----------------------------------------------------------------------------
vtkOGSSpatialStats::~vtkOGSSpatialStats() {
	this->StatDataArraySelection->Delete();

	#ifdef PARAVIEW_USE_MPI
		this->SetController(NULL);	
	#endif
}

//----------------------------------------------------------------------------
int vtkOGSSpatialStats::RequestInformation(vtkInformation* vtkNotUsed(request),
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

//----------------------------------------------------------------------------
int vtkOGSSpatialStats::RequestData(vtkInformation *vtkNotUsed(request),
	vtkInformationVector **inputVector,vtkInformationVector *outputVector) {
	// Get the info objects
	vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
	vtkInformation *outInfo = outputVector->GetInformationObject(0);

	// Stop all threads except from the master to execute
	#ifdef PARAVIEW_USE_MPI
	if (this->procId > 0) return 1;
	#endif

	// Get the input and output
	vtkDataSet *input = vtkDataSet::SafeDownCast(
		inInfo->Get(vtkDataObject::DATA_OBJECT()));
	vtkDataSet *output = vtkDataSet::SafeDownCast(
		outInfo->Get(vtkDataObject::DATA_OBJECT()));

	// We just want to copy the mesh, not the variables
	output->CopyStructure(input);
	output->GetFieldData()->PassData(input->GetFieldData());

	this->UpdateProgress(0.);

	// Test for cell data or point data
	if (VTKARRAY::SafeDownCast(input->GetCellData()->GetArray("e1")) == NULL) {
		this->iscelld  = false; 
	}

	// This section is only executed once, to populate the xyz and
	// cId2zId arrays. Successive iterations should not execute.
	// This section is included here since RequestInformation gives
	// troubles when restarting.
	if (this->xyz.isempty() || this->isReqInfo) {

		this->isReqInfo = false;
		
		// Recover Metadata array (depth factor)
		vtkStringArray *vtkmetadata = vtkStringArray::SafeDownCast(
			input->GetFieldData()->GetAbstractArray("Metadata"));
		double dfact = (vtkmetadata != NULL) ? std::stod( vtkmetadata->GetValue(2) ) : 1000.;
		
		if (vtkmetadata == NULL) 
			vtkWarningMacro("Field array Metadata not found! Depth factor set to 1000. automatically.");

		// Recover cell or point coordinates, depending on the array type
		this->xyz = (this->iscelld) ? VTK::getVTKCellCenters(input,dfact) : 
		                              VTK::getVTKCellPoints(input,dfact);

		// If the user has not inputed any depth level, make sure to clear the vector
		if(this->ndepths == 0) this->zcoords.clear();

		// Up to this point we have the cell centers or point coordinates correctly
		// stored under "xyz". Now we shall find the number of unique z coordinates or,
		// depending on the user input, the coordinates of each depth level, as well as
		// its mesh connectivity (cId2zId).
		this->cId2zId = field::countDepthLevels(this->xyz,this->zcoords,this->epsi,true);
	}

	// At this point we either have a rectilinear grid or an
	// unstructured grid with either cell or point data
	//
	// Use "e1" to try to discern whether we have cell data or
	// point data. Since "e1", "e2" and "e3" are needed, we cannot
	// compute unless they are present as either cell or point data
	VTKARRAY *vtke1 = NULL, *vtke2 = NULL, *vtke3 = NULL;

	vtke1 = (this->iscelld) ? VTKARRAY::SafeDownCast(input->GetCellData()->GetArray("e1")) :
						      VTKARRAY::SafeDownCast(input->GetPointData()->GetArray("e1"));
	vtke2 = (this->iscelld) ? VTKARRAY::SafeDownCast(input->GetCellData()->GetArray("e2")) :
						      VTKARRAY::SafeDownCast(input->GetPointData()->GetArray("e2"));
	vtke3 = (this->iscelld) ? VTKARRAY::SafeDownCast(input->GetCellData()->GetArray("e3")) :
						      VTKARRAY::SafeDownCast(input->GetPointData()->GetArray("e3"));

	if (vtke1 == NULL || vtke2 == NULL || vtke3 == NULL) {
		vtkErrorMacro("Mesh weights (e1, e2 and e3) need to be loaded to proceed!");
		return 0;
	}

	// Convert to field arrays
	field::Field<FLDARRAY> e1, e2, e3;
	e1 = VTK::createFieldfromVTK<VTKARRAY,FLDARRAY>(vtke1);
	e2 = VTK::createFieldfromVTK<VTKARRAY,FLDARRAY>(vtke2);
	e3 = VTK::createFieldfromVTK<VTKARRAY,FLDARRAY>(vtke3);

	this->UpdateProgress(.1);

	// Output connectivity (Uncomment for debugging purposes)
	//vtkIntArray *vtkcId2zId;
	//vtkcId2zId = VTK::createVTKfromField<vtkIntArray,int>("cId2zId",this->cId2zId);
	//output->GetCellData()->AddArray(vtkcId2zId);
	//vtkcId2zId->Delete();

	// We can now loop the number of active variables
	int narrays = (this->iscelld) ? input->GetCellData()->GetNumberOfArrays() :
									input->GetPointData()->GetNumberOfArrays();
	int nstat   = this->GetNumberOfStatArrays();

	// Parallelization strategy MPI
	for (int varId = 0; varId < narrays; ++varId) {
		// Recover the array and the array name
		vtkDataArray *vtkDArray;
		vtkDArray = (this->iscelld) ? input->GetCellData()->GetArray(varId) :
									  input->GetPointData()->GetArray(varId);
		std::string arrName = vtkDArray->GetName();

		// We should not average the coast or basins mask nor e1t or e2t
		// Names have been harcoded here as there is no way to ensure that
		// these arrays will exist or not.
		if (arrName == std::string("coast mask"))  continue;
		if (arrName == std::string("basins mask")) continue;
		if (arrName == std::string("land mask"))   continue;
		if (arrName == std::string("e1"))          continue;
		if (arrName == std::string("e2"))          continue;
		if (arrName == std::string("e3"))          continue;

		// Recover Array values
		VTKARRAY *vtkArray;
		field::Field<FLDARRAY> array;
		vtkArray = VTKARRAY::SafeDownCast( vtkDArray );
		array = VTK::createFieldfromVTK<VTKARRAY,FLDARRAY>(vtkArray);

		// Also anything not being an scalar array should not be computed
		if (array.get_m() > 1) {
			vtkWarningMacro("Variable " << arrName.c_str() <<
				" is not an scalar array. Statistics will be skipped for this variable.");
			continue;
		}

		// We need to decide which statistics to compute. Create a map associating
		// the statistic with a field. Initialize to zero.
		std::map<std::string,field::Field<FLDARRAY>> mStatArray;
		for (int statId = 0; statId < nstat; ++statId) {
			// Recover stat name
			std::string statName = this->GetStatArrayName(statId);
			// Skip those arrays who have not been enabled
			if (!this->GetStatArrayStatus(statName.c_str()))
				continue;
			// Create empty field array
			field::Field<FLDARRAY> statArray(array.get_n(),array.get_m(),0.);
			// Pair and insert to map
			mStatArray.insert(std::make_pair(statName,statArray));
		}

		// Also, as a preprocess, create some maps to refer the data on the layers
		// to the data on the mesh. Populate them by looping the mesh once.
		std::map<int,std::vector<FLDARRAY>> mValuesPerLayer, mWeightPerLayer;
		std::map<int,std::vector<int>> mMeshPerLayer;

		#pragma omp parallel
		{
		for (int ii=OMP_THREAD_NUM; ii < array.get_n(); ii+=OMP_NUM_THREADS) {
			// Set maps
			double w = (this->useVolume) ? e1[ii][0]*e2[ii][0]*e3[ii][0] : e1[ii][0]*e2[ii][0];
			int ind  = cId2zId[ii][0];
			#pragma omp critical
			{
			mValuesPerLayer[ind].push_back(array[ii][0]);
			mWeightPerLayer[ind].push_back(w);
			mMeshPerLayer[ind].push_back(ii);
			}	
		}
		
		// Wait until the previous calculation has been performed
		// then proceed to run for each depth layer
		#pragma omp barrier
		
		// For each depth layer, compute the statistics
		for (int zId = OMP_THREAD_NUM; zId < this->zcoords.size(); zId += OMP_NUM_THREADS) {

			std::vector<std::pair<FLDARRAY,int>> sortedValues(mValuesPerLayer[zId].size());

			// Iterate on the layer, compute  the weights, mean, std dev and min/max
			// This code uses the weighted variant of Welford's online algorithm for
			// the standard deviation.
			// West, D. H. D. "Updating mean and variance estimates: An improved method." Communications of the ACM 22.9 (1979): 532-535.
			double sum_weight = 0., meanval = 0., meanval_old = 0., stdval = 0., maxval = -1.e20, minval = 1.e20;

			for (int ii=0; ii<mValuesPerLayer[zId].size(); ++ii) {
				double v = mValuesPerLayer[zId][ii];
				double w = mWeightPerLayer[zId][ii];
				// Minimum and maximum
				minval = (v < minval) ? v : minval;
				maxval = (v > maxval) ? v : maxval;
				// Weights, mean and standard deviation
				sum_weight += w;
				meanval_old = meanval;
				meanval    += w/sum_weight*(v-meanval);
				stdval     += w*(v-meanval_old)*(v-meanval);
				// Array to sort
				sortedValues[ii] = std::make_pair(v,ii);				
			}

			// Finish computing standard deviation
			stdval = sqrt(stdval/sum_weight);

			// Percentile Weight computation
			// Just execute this part of the code if the user has requested to compute the
			// percentiles. Otherwise, save a couple of for loops.
			double perc[]    = {.05,.25,.50,.75,.95};
			double percval[] = { 0., 0., 0., 0., 0.};	
			
			if ( this->GetStatArrayStatus("p05") || this->GetStatArrayStatus("p50") ||
				 this->GetStatArrayStatus("p75") || this->GetStatArrayStatus("p95") ) {

				// Order the values
				std::sort(sortedValues.begin(),sortedValues.end());

				// Second iteration on the layer, this time compute the percentile weight
				double pweight = 0.;
				std::vector<FLDARRAY> percw(mValuesPerLayer[zId].size());
				
				for (int ii=0; ii<mValuesPerLayer[zId].size(); ++ii) {
					int ind  = sortedValues[ii].second;
					pweight += mWeightPerLayer[zId][ind];
					double aux = (pweight - .5*mWeightPerLayer[zId][ind])/sum_weight; // Reused variable
					percw[ii] = aux;
				}

				// Compute the percentiles
				for (int pp = 0; pp < 5; pp++) {
					// Find the value that is equal to perc or immediately after.
					std::vector<FLDARRAY>::iterator lbound = std::lower_bound(percw.begin(),percw.end(),perc[pp]);
					// This is our position on the ordered value array
					int s = (lbound - percw.begin()) < 0 ? 0 : lbound - percw.begin(); 
					// Set the value for the weight
					if (s == 0)                             { percval[pp] = sortedValues[s].first; continue; } // == sd[0]
					if (s == mValuesPerLayer[zId].size()-1) { percval[pp] = sortedValues[s].first; continue; } // == sd[n-1]
					double f1 = (percw[s] - perc[pp])   / (percw[s] - percw[s-1]);
					double f2 = (perc[pp] - percw[s-1]) / (percw[s] - percw[s-1]);

					percval[pp] = f1*sortedValues[s-1].first + f2*sortedValues[s].first;
				}
			}

			// Third iteration on the layers, this time set
			// the mesh to the proper value.
			for (int ii=0; ii<mMeshPerLayer[zId].size(); ++ii) {
				std::map<std::string,field::Field<FLDARRAY>>::iterator iter;
				int ind = mMeshPerLayer[zId][ii];
				// Loop each one of the active statistics
				for (iter = mStatArray.begin(); iter != mStatArray.end(); ++iter) {
					// Switch computed statistic
					switch ( this->GetStatArrayIndex(iter->first.c_str()) ) {
						case 0: // Mean
							iter->second[ind][0] = meanval; break;
						case 1: // Std dev
							iter->second[ind][0] = stdval; break;
						case 2: // Min
							iter->second[ind][0] = minval; break;
						case 3: // p05
							iter->second[ind][0] = percval[0]; break;;
						case 4: // p25
							iter->second[ind][0] = percval[1]; break;
						case 5: // p50
							iter->second[ind][0] = percval[2]; break;
						case 6: // p75
							iter->second[ind][0] = percval[3]; break;
						case 7: // p95
							iter->second[ind][0] = percval[4]; break;
						case 8: // Max
							iter->second[ind][0] = maxval; break;
					}
				}				
			}
		}
		}

		// Now that we computed the arrays, we can set them in the output
		// and deallocate memory
		std::map<std::string,field::Field<FLDARRAY>>::iterator iter;
		for (iter = mStatArray.begin(); iter != mStatArray.end(); ++iter) {
			// Generate variable name
			std::string statVarName = arrName + std::string(", ") + iter->first;
			// Convert field to VTK
			vtkArray = VTK::createVTKfromField<VTKARRAY,FLDARRAY>(statVarName.c_str(),iter->second);
			if (this->iscelld)
				output->GetCellData()->AddArray(vtkArray);
			else
				output->GetPointData()->AddArray(vtkArray);
			vtkArray->Delete();
		}

		this->UpdateProgress(0.1+0.9/(double)(narrays)*(double)(varId));
	}

	this->UpdateProgress(1.);
	return 1;
}

//----------------------------------------------------------------------------
void vtkOGSSpatialStats::SetNumberOfDepthLevels(int n)
{
	// We basically use this function to set nz
	this->ndepths = n > 0 ? n+1 : 0;
	// Also restart the vector that stores the coordinates
	this->zcoords.clear();
	this->Modified();
}

void vtkOGSSpatialStats::SetDepthLevels(int i, double value)
{
	// We basically use this function to set zcoords
	this->zcoords.push_back(-value);
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
