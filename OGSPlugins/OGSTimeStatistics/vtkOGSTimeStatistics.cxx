/*=========================================================================

  Program:   OGSTimeStatistics
  Module:    vtkOGSTimeStatistics.cxx

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkOGSTimeStatistics.h"

#include "vtkDataSet.h"
#include "vtkCellData.h"
#include "vtkPointData.h"
#include "vtkFieldData.h"
#include "vtkStringArray.h"
#include "vtkFloatArray.h"
#include "vtkDoubleArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#include "vtkObjectFactory.h"

#include <cmath>
#include <ctime>
#include <chrono>
#include <vector>
#include <string>

#ifdef __GNUC__
// Include OpenMP when working with GCC
#include <omp.h>
#define OMP_NUM_THREADS omp_get_num_threads()
#define OMP_THREAD_NUM  omp_get_thread_num()
#else
#define OMP_NUM_THREADS 1
#define OMP_THREAD_NUM  0
#endif

vtkStandardNewMacro(vtkOGSTimeStatistics);

#ifdef PARAVIEW_USE_MPI
#include "vtkMultiProcessController.h"
vtkCxxSetObjectMacro(vtkOGSTimeStatistics, Controller, vtkMultiProcessController);
#endif

//----------------------------------------------------------------------------

/*
	Macro to set the array precision 
*/
#define FLDARRAY double
#define VTKARRAY vtkDoubleArray

#include "../_utils/field.h"
#include "../_utils/vtkFields.hpp"
#include "../_utils/OGS.hpp"
#include "../_utils/netcdfio.hpp"

//----------------------------------------------------------------------------
void InitializeStatistics(std::vector<std::string> &arrNames, std::vector<std::string> &arrNamesNotComputed,
	vtkRectilinearGrid *input, vtkRectilinearGrid *output);
void RecoverMasterFileName(std::string &fname, vtkRectilinearGrid *input);

//----------------------------------------------------------------------------
vtkOGSTimeStatistics::vtkOGSTimeStatistics() {
	this->TimeValues = vtkStringArray::New();

	this->ii_start = 0;
	this->ii_end   = 0;

	this->procId = 0;
	this->nProcs = 0;

	#ifdef PARAVIEW_USE_MPI
	this->Controller = nullptr;
	this->SetController(vtkMultiProcessController::GetGlobalController());
	#endif
}

//----------------------------------------------------------------------------
vtkOGSTimeStatistics::~vtkOGSTimeStatistics() {
	this->TimeValues->Delete();

	#ifdef PARAVIEW_USE_MPI
	this->SetController(nullptr);
	#endif
}

//----------------------------------------------------------------------------
int vtkOGSTimeStatistics::RequestInformation(vtkInformation* vtkNotUsed(request),
  vtkInformationVector** vtkNotUsed(inputVector), vtkInformationVector* outputVector) {

	// Grab the input
	vtkInformation *outInfo = outputVector->GetInformationObject(0);

	// Handle the parallel controller
	#ifdef PARAVIEW_USE_MPI
	if (this->Controller->GetNumberOfProcesses() > 1) {
		this->procId = this->Controller->GetLocalProcessId();
		this->nProcs = this->Controller->GetNumberOfProcesses();
	}
	#endif

	// Populate the TimeValues array. For that we will use the data
	// stored in the timesteps to be consistent and since the metadata
	// array might not be available from the beginning.
	int ntsteps = outInfo->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
	double *tsteps; tsteps = outInfo->Get(vtkStreamingDemandDrivenPipeline::TIME_STEPS());

	char buff[256];
	this->TimeValues->Delete();
	this->TimeValues = VTK::createVTKstrf("TimeValues",ntsteps,NULL);
	
	for (int ii = 0; ii < ntsteps; ++ii) {
		// Convert to struct tm
		time_t time = std::chrono::system_clock::to_time_t(std::chrono::system_clock::time_point(
			std::chrono::duration_cast<std::chrono::seconds>(std::chrono::duration<double>(tsteps[ii]))));
		struct tm tm = *localtime(&time);

		// Format and display
		strftime(buff,256,"%Y%m%d-%H:%M:%S",&tm);
		this->TimeValues->SetValue(ii,buff);
	}

	// The output data of this filter has no time associated with it. It is the
	// result of computations that happen over all time.
	outInfo->Remove(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
	outInfo->Remove(vtkStreamingDemandDrivenPipeline::TIME_RANGE());

	return 1;
}

//----------------------------------------------------------------------------
int vtkOGSTimeStatistics::RequestData(vtkInformation *request, 
	vtkInformationVector **inputVector, vtkInformationVector *outputVector) {

	// get the info objects
	vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
	vtkInformation *outInfo = outputVector->GetInformationObject(0);

	// get the input and output
	vtkRectilinearGrid *input = vtkRectilinearGrid::SafeDownCast(
		inInfo->Get(vtkDataObject::DATA_OBJECT()));
	vtkRectilinearGrid *output = vtkRectilinearGrid::SafeDownCast(
		outInfo->Get(vtkDataObject::DATA_OBJECT()));

	// Copy the mesh, not the variables
	output->CopyStructure(input);

	/* INITIALIZATION PHASE

		Obtain the names of the variables that need to be time-averaged by reading
		which variables have been loaded as cell arrays. Note that this filter cannot
		deal with point array variables.

		Recover (and broadcast) the file name of the OGS master file.

		Initialize arrays for accumulating.

	*/
	// Define a vector containing the names of the arrays to compute and to exclude
	std::string FileName;
	std::vector<std::string> arrNames, arrNamesNotComputed{"e1","e2","e3",
		"basins mask","coast mask","Okubo-Weiss mask","Q-criterion mask"};

	VTKARRAY *vtkArray;
	std::vector<field::Field<FLDARRAY>> arrFields;
	
	#ifdef PARAVIEW_USE_MPI
	
	// Thread 0 has all the information in input and output, therefore is the
	// only computing the array which will later be broadcasted to all ranks
	if (this->procId == 0) {
		InitializeStatistics(arrNames,arrNamesNotComputed,input,output);
		RecoverMasterFileName(FileName, input);
		// Initalize accumulating fields
		for (std::string arrName : arrNames) {
			vtkArray = VTKARRAY::SafeDownCast( input->GetCellData()->GetArray(arrName.c_str()) );
			arrFields.push_back( field::Field<FLDARRAY>(vtkArray->GetNumberOfTuples(),vtkArray->GetNumberOfComponents(),0.) );
		}
	}
	
	// Broadcast the information to all other threads if applicable
	if (this->nProcs > 1) {
		// Broadcast master file name
		int str_len = FileName.length();
		this->Controller->Broadcast(&str_len,1,0);
		char buff[128] = ""; sprintf(buff,"%s",FileName.c_str());
		this->Controller->Broadcast(buff,str_len,0);
		FileName = std::string(buff);
		
		// Broadcast number of arrays to compute
		int n_arrays = arrNames.size();
		this->Controller->Broadcast(&n_arrays,1,0);
		
		// Allocate array on the other ranks
		if (this->procId > 0) arrNames.resize(n_arrays,std::string(""));
		for (int varId = 0; varId < arrNames.size(); ++varId) {
			// Array length
			str_len = arrNames[varId].length();
			this->Controller->Broadcast(&str_len,1,0);
			// Array name
			std::memset(buff,0,128*sizeof(char));
			sprintf(buff,"%s",arrNames[varId].c_str());
			this->Controller->Broadcast(buff,str_len,0);
			arrNames[varId] = std::string(buff);
			// Array field
			int n = 0, m = 0;
			vtkArray = VTKARRAY::SafeDownCast( input->GetCellData()->GetArray(arrNames[varId].c_str()) );
			if (vtkArray) {
				n = vtkArray->GetNumberOfTuples();
				m = vtkArray->GetNumberOfComponents();
			}
			this->Controller->Broadcast(&n,1,0);
			this->Controller->Broadcast(&m,1,0);
			if (this->procId > 0) 
				arrFields.push_back( field::Field<FLDARRAY>(n,m,0.) );
		}
	}
	
	#else
	
	// This is the normal non-parallel algorithm
	InitializeStatistics(arrNames,arrNamesNotComputed,input,output);
	RecoverMasterFileName(FileName, input);

	// Initalize accumulating fields
	for (std::string arrName : arrNames) {
		vtkArray = VTKARRAY::SafeDownCast( input->GetCellData()->GetArray(arrName.c_str()) );
		arrFields.push_back( field::Field<FLDARRAY>(vtkArray->GetNumberOfTuples(),vtkArray->GetNumberOfComponents(),0.) );
	}

	#endif

	// Read the OGS file to be able to generate the paths to the variable files
	ogs::OGS ogsdata(FileName);
	if (ogsdata.readMainFile() == -1) { 
		vtkErrorMacro("Cannot read <"<<FileName<<">!\nAborting");
		return 0;
	}

	this->UpdateProgress(0.1);

	/* ACCUMULATE PHASE

		Open and read the NetCDF files containing the variables that we previously
		defined. Accumulate in the arrays from ii_start to ii_end.

	*/
	
	int time_interval[2] = {this->ii_start,this->ii_end};

	// Parallel partition
	#ifdef PARAVIEW_USE_MPI

	if (this->nProcs > 1) {
		// Main thread (0) contains values for ii_start and ii_end
		// They must be sent to the other processes
		this->Controller->Broadcast(time_interval,2,0);

		// Now everyone should have the range of start and end times
		// We must split equally among the threads. We must also control
		// that the number of threads is less or equal than the intervals
		// requested.
		int range = time_interval[1] - time_interval[0];
		if (this->nProcs < range) {
			// We split normally among processes assuming no remainder
			int rangePerProcess = std::floor(range/this->nProcs);
			this->ii_start = time_interval[0] + this->procId*rangePerProcess;
			this->ii_end   = this->ii_start + rangePerProcess;
			// Handle the remainder
			int remainder = range - rangePerProcess*this->nProcs;
			if (remainder > this->procId){
				this->ii_start += this->procId;
				this->ii_end   += this->procId + 1;
			} else {
				this->ii_start += remainder;
				this->ii_end   += remainder;
			}
		} else {
			// Each process will forcefully conduct one instant.
			this->ii_start = (this->procId < time_interval[1]) ? this->procId     : time_interval[1];
			this->ii_end   = (this->procId < time_interval[1]) ? this->procId + 1 : time_interval[1];
		}
	}

	#endif

	// Ensure the validity of the time range
	if (time_interval[0] > time_interval[1]) {
		if (this->procId == 0) 
			vtkErrorMacro("End time is greater than initial time! Please select a valid time range.");
		return 0;
	}

	// Loop the instants, prepare arrays and variables for the computation of the mean
	field::Field<FLDARRAY> arrayTemp;
	
	// For each timestep
	for(int ii = this->ii_start; ii < this->ii_end; ++ii) {
		FLDARRAY ii_range = (FLDARRAY)(ii + 1. - this->ii_start);
		// For each variable 
		for (int varId = 0; varId < arrNames.size(); ++varId) {
			if (arrNames[varId] == std::string("Velocity")) {
				// Load the variable on a temporal field
				arrayTemp.set_dim(arrFields[varId].get_n(),arrFields[varId].get_m());

				if ( NetCDF::readNetCDF2F3(ogsdata.var_path(arrNames[varId],ii).c_str(),
											"vozocrtx","vomecrty","vovecrtz", arrayTemp ) != NETCDF_OK )	{					    								  
					vtkErrorMacro("Cannot read variable <"<<arrNames[varId]<<"> in NetCDF! Aborting!"); return 0;
				}
				// Projecting the velocity field to the UVW grid is done a the end
			} else {
				// Load the variable on a temporal field
				arrayTemp.set_dim(arrFields[varId].get_n(),arrFields[varId].get_m());

				if ( NetCDF::readNetCDF2F(ogsdata.var_path(arrNames[varId],ii).c_str(), 
					ogsdata.var_vname(arrNames[varId]), arrayTemp) != NETCDF_OK ) {
					vtkErrorMacro("Cannot read variable <"<<arrNames[varId]<<"> in NetCDF! Aborting!"); return 0;
				}
			}
			// Accumulate, implement 1 pass averaging algorithm
			#pragma omp parallel shared(arrFields,arrayTemp) firstprivate(ii_range)
			{
			for (int nId = OMP_THREAD_NUM; nId < arrFields[varId].get_n(); nId += OMP_NUM_THREADS) {
				for (int mId = 0; mId < arrFields[varId].get_m(); ++mId) { 
					arrFields[varId][nId][mId] += (1./ii_range)*(arrayTemp[nId][mId] - arrFields[varId][nId][mId]);
				}
			}
			}
			arrayTemp.clear();
		}
		this->UpdateProgress(0.1+0.7/(this->ii_end-this->ii_start)*(ii-this->ii_start));
	}

	/* REDUCTION PHASE

		If run in parallel with more than one rank, 
		the reduction of the fields is done here.

	*/
	#ifdef PARAVIEW_USE_MPI

	// Only reduce if we have more than 1 process
	if (this->nProcs > 1) {
		FLDARRAY n_range = (FLDARRAY)(this->ii_end-this->ii_start);
		// For each variable 
		for (int varId = 0; varId < arrNames.size(); ++varId) {
			arrayTemp.set_dim(arrFields[varId].get_n(),arrFields[varId].get_m());

			// Reduce to arrayTemp on rank 0
			arrFields[varId] *= n_range;
			this->Controller->Reduce(arrFields[varId].data(),arrayTemp.data(),arrayTemp.get_sz(),
				vtkCommunicator::StandardOperations::SUM_OP,0);

			// Average among processes
			if (this->procId == 0) {
				#pragma omp parallel shared(arrFields,arrayTemp) firstprivate(time_interval)
				{
				for (int nId = OMP_THREAD_NUM; nId < arrFields[varId].get_n(); nId += OMP_NUM_THREADS) {
					for (int mId = 0; mId < arrayTemp.get_m(); ++mId) { 
						arrFields[varId][nId][mId] = arrayTemp[nId][mId]/( (FLDARRAY)(time_interval[1] - time_interval[0]) );
					}
				}
				}
			}

			arrayTemp.clear();
		}
	}

	#endif

	this->UpdateProgress(0.9);

	/* FINALIZE

		Finalization, the master process stores the arrays
		inside the output.

	*/

	// Add arrays to output
	if (this->procId == 0) {
		for (int varId = 0; varId < arrNames.size(); ++varId) {
			vtkArray = VTK::createVTKfromField<VTKARRAY,FLDARRAY>(arrNames[varId],arrFields[varId]);
			output->GetCellData()->AddArray(vtkArray);
			vtkArray->Delete();
		}
		// Make sure the metadata array is passed to the output
		vtkStringArray *vtkmetadata = vtkStringArray::SafeDownCast(
			input->GetFieldData()->GetAbstractArray("Metadata"));
		output->GetFieldData()->AddArray(vtkmetadata);
	}

	this->UpdateProgress(1.);
	return 1;
}

//----------------------------------------------------------------------------
void InitializeStatistics(std::vector<std::string> &arrNames, std::vector<std::string> &arrNamesNotComputed,
	vtkRectilinearGrid *input, vtkRectilinearGrid *output) {

	// We will now copy cell  arrays to the arrNames, excluding
	// these on arrNamesNotComputed. These not computed will be
	// added to the output

	vtkDataArray *vtkDArray;
	arrNames.clear();

	for (int varId = 0; varId < input->GetCellData()->GetNumberOfArrays(); ++varId) {
		// Recover the array and the array name
		vtkDArray = input->GetCellData()->GetArray(varId);
		std::string arrName = vtkDArray->GetName();

		// Exclude arrNamesNotComputed and add them to output
		bool addArray = true;
		for (std::string arrNameNotComputed : arrNamesNotComputed) {
			if (arrNameNotComputed == arrName) {
				output->GetCellData()->AddArray(vtkDArray);
				addArray = false;
				break;
			}
		}

		if (addArray) arrNames.push_back(arrName);
	}
}

//----------------------------------------------------------------------------
void RecoverMasterFileName(std::string &fname, vtkRectilinearGrid *input) {

	// Recover the master file name from the metadata array
	// Return whether we need to stop executing or not

	vtkStringArray *vtkmetadata = vtkStringArray::SafeDownCast(
		input->GetFieldData()->GetAbstractArray("Metadata"));

	// If successful, recover the file name
	if (vtkmetadata)
		fname = vtkmetadata->GetValue(4);
}

//----------------------------------------------------------------------------
void vtkOGSTimeStatistics::SetStartTime(const char *tstep) {
	// Obtain the timestep index
	for (int ii = 0; ii < this->TimeValues->GetNumberOfTuples(); ++ii) {
		if (std::string(this->TimeValues->GetValue(ii)) == std::string(tstep)) {
			this->ii_start = ii; break;
		}
	}
	this->Modified();
}

void vtkOGSTimeStatistics::SetEndTime(const char *tstep) {
	// Obtain the timestep index
	for (int ii = 0; ii < this->TimeValues->GetNumberOfTuples(); ++ii) {
		if (std::string(this->TimeValues->GetValue(ii)) == std::string(tstep)) {
			this->ii_end = ii; break;
		}
	}
	this->Modified();
}

//----------------------------------------------------------------------------
vtkStringArray *vtkOGSTimeStatistics::GetTimeValues() {
	return this->TimeValues;
}
