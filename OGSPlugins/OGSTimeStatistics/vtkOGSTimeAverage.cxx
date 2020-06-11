/*=========================================================================

  Program:   OGSTimeAverage
  Module:    vtkOGSTimeAverage.cxx

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

	 This software is distributed WITHOUT ANY WARRANTY; without even
	 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
	 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkOGSTimeAverage.h"

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

vtkStandardNewMacro(vtkOGSTimeAverage);

#ifdef PARAVIEW_USE_MPI
#include "vtkMultiProcessController.h"
vtkCxxSetObjectMacro(vtkOGSTimeAverage, Controller, vtkMultiProcessController);
#endif

//----------------------------------------------------------------------------
#include "macros.h"
#include "field.h"
#include "TimeObject.h"
#include "TimeInterval.h"
#include "TimeList.h"
//#include "netcdfio.h"
//#include "OGS.hpp"
#include "vtkFields.h"

//----------------------------------------------------------------------------
void BuildTimeList(Time::TimeList &TL, vtkInformation *Info) {

	// Build a TimeList object using the pipeline temporal data.
	// This TimeList will be later used for computing the averages
	// given a TimeRequestor. The metadata array might not be available 
	// from the beginning.

	int ntsteps = Info->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
	double *tsteps; tsteps = Info->Get(vtkStreamingDemandDrivenPipeline::TIME_STEPS());

	// Create a TimeObjectList where to store all the instants
	Time::TimeObjectList TOL(ntsteps);

	// Iterate the number of steps and set the values of the list
	for (int ii = 0; ii < ntsteps; ++ii) {
		// Convert to struct tm
		time_t time = std::chrono::system_clock::to_time_t(std::chrono::system_clock::time_point(
			std::chrono::duration_cast<std::chrono::seconds>(std::chrono::duration<double>(tsteps[ii]))));
		struct tm tm = *localtime(&time);

		// Set up the TimeObjectList
		TOL[ii] = Time::TimeObject(tm);
	}

	// Sort and print the list (for debugging purposes)
	TOL.sort();
//	printf("TOL: %s\n",TOL.as_string("%Y%m%d").c_str());

	// Now create the TimeList from the TimeObjectList
	TL = Time::TimeList(TOL);
//	printf("Defined list of %d elements: %s ... %s\n",TL.len(),TL[0].as_string("%Y-%m-%d %H:%M:%S").c_str(),TL[-1].as_string("%Y-%m-%d %H:%M:%S").c_str());
}

////----------------------------------------------------------------------------
//void InitializeStatistics(std::vector<std::string> &arrNames, std::vector<std::string> &arrNamesNotComputed,
//	vtkRectilinearGrid *input, vtkRectilinearGrid *output) {
//
//	// We will now copy cell  arrays to the arrNames, excluding
//	// these on arrNamesNotComputed. These not computed will be
//	// added to the output
//
//	vtkDataArray *vtkDArray;
//	arrNames.clear();
//
//	for (int varId = 0; varId < input->GetCellData()->GetNumberOfArrays(); ++varId) {
//		// Recover the array and the array name
//		vtkDArray = input->GetCellData()->GetArray(varId);
//		std::string arrName = vtkDArray->GetName();
//
//		// Exclude arrNamesNotComputed and add them to output
//		bool addArray = true;
//		for (std::string arrNameNotComputed : arrNamesNotComputed) {
//			if (arrNameNotComputed == arrName) {
//				output->GetCellData()->AddArray(vtkDArray);
//				addArray = false;
//				break;
//			}
//		}
//
//		if (addArray) arrNames.push_back(arrName);
//	}
//}
//
////----------------------------------------------------------------------------
//void RecoverMasterFileName(std::string &fname, vtkRectilinearGrid *input) {
//
//	// Recover the master file name from the metadata array
//	// Return whether we need to stop executing or not
//
//	vtkStringArray *vtkmetadata = vtkStringArray::SafeDownCast(
//		input->GetFieldData()->GetAbstractArray("Metadata"));
//
//	// If successful, recover the file name
//	if (vtkmetadata)
//		fname = vtkmetadata->GetValue(4);
//}

//----------------------------------------------------------------------------
vtkOGSTimeAverage::vtkOGSTimeAverage() {
	this->sum_weights      = 0.;
	this->TL_computed      = false;
	this->use_files        = false;

	this->CurrentTimeIndex = 0;

	this->procId = 0;
	this->nProcs = 0;

	#ifdef PARAVIEW_USE_MPI
	this->Controller = nullptr;
	this->SetController(vtkMultiProcessController::GetGlobalController());
	#endif
}

//----------------------------------------------------------------------------
vtkOGSTimeAverage::~vtkOGSTimeAverage() {

	#ifdef PARAVIEW_USE_MPI
	this->SetController(nullptr);
	#endif
}

//----------------------------------------------------------------------------
int vtkOGSTimeAverage::RequestInformation(vtkInformation* vtkNotUsed(request),
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

	// Compute the TimeList if it hasn't been computed
	if (!this->TL_computed) {
		BuildTimeList(this->TL,outInfo);
		this->TL_computed = true;
	}

	// The output data of this filter has no time associated with it. It is the
	// result of computations that happen over all time.
	outInfo->Remove(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
	outInfo->Remove(vtkStreamingDemandDrivenPipeline::TIME_RANGE());

	return 1;
}

//------------------------------------------------------------------------------
int vtkOGSTimeAverage::RequestUpdateExtent(vtkInformation* vtkNotUsed(request),
	vtkInformationVector** inputVector, vtkInformationVector* vtkNotUsed(outputVector)) {
	
	vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);

	// The RequestData method will tell the pipeline executive to iterate the
	// upstream pipeline to get each time step in order.  The executive in turn
	// will call this method to get the extent request for each iteration (in this
	// case the time step).
	double* inTimes = inInfo->Get(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
	
	if (inTimes) {
		int inst = (this->instants.size() == 0) ? this->CurrentTimeIndex : this->instants[this->CurrentTimeIndex];
//		printf("Time: <%d> %s\n",inst,this->TL[inst].as_string("%Y-%m-%d %H:%M:%S").c_str());
		inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP(), inTimes[inst]);
	}
	
	return 1;
}

//----------------------------------------------------------------------------
int vtkOGSTimeAverage::RequestData(vtkInformation *request, 
	vtkInformationVector **inputVector, vtkInformationVector *outputVector) {

	// get the info objects
	vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
	vtkInformation *outInfo = outputVector->GetInformationObject(0);

	// get the input and output
	vtkDataSet *input = vtkDataSet::SafeDownCast(
		inInfo->Get(vtkDataObject::DATA_OBJECT()));
	vtkDataSet *output = vtkDataSet::SafeDownCast(
		outInfo->Get(vtkDataObject::DATA_OBJECT()));

	if (this->CurrentTimeIndex == 0) {
		// Master computes the instants to loop
		if (this->procId == 0) {
			// Compute the TimeList if it hasn't been computed
			if (!this->TL_computed) {
				BuildTimeList(this->TL,outInfo);
				this->TL_computed = true;
			}
//			printf("Defined list of %d elements: %s ... %s\n",this->TL.len(),
//				this->TL[0].as_string("%Y-%m-%d %H:%M:%S").c_str(),
//				this->TL[-1].as_string("%Y-%m-%d %H:%M:%S").c_str());

			// At this point check that the TimeInterval length is not zero
			// if so, crash and stop, there is no sense to compute an average.
//			printf("Time interval is %s of length %ld\n",
//				this->TI.as_string("%Y%m%d-%H:%M:%S").c_str(),this->TI.length());
			if (this->TI.length() == 0) {
				vtkErrorMacro("TimeInterval is <"<<this->TI.as_string("%Y%m%d-%H:%M:%S")<<">! Cannot continue...");
				return 0;
			}

			// We can now create a GenericRequestor with this TimeInterval.
			// We will use it to obtain all the instants and average weights.
			Time::Requestor req(this->TI,"%Y%m%d-%H:%M:%S");

			if (this->TL.select(&req,this->instants,this->weights) == TIME_ERR) {
				vtkErrorMacro("Error selecting instants with Requestor! Cannot continue...");
				return 0;
			}
//			for (int ii=0; ii<instants.size();++ii)
//				printf("ind: %d, weight: %f, %s\n",this->instants[ii],this->weights[ii],
//					this->TL[instants[ii]].as_string("%Y-%m-%d %H:%M:%S").c_str());
		}

		// Broadcast instants and weights to all the processors if needed
		if (this->nProcs > 1) {
			// Array length
			int arr_len = (int)(this->instants.size());
			this->Controller->Broadcast(&arr_len,1,0);
			// Allocate
			if (this->procId > 0) { this->instants.resize(arr_len,0); this->weights.resize(arr_len,0.); }
			// Broadcast array
			this->Controller->Broadcast(this->instants.data(),arr_len,0);
			this->Controller->Broadcast(this->weights.data(),arr_len,0);
		}

		// Compute the sum of weights for the average
		this->sum_weights = 0.;
		for (int ii=0; ii<instants.size();++ii)
			this->sum_weights += this->weights[ii];
	}

	if (this->use_files) { 
		// Optimized algorithm based on parallel access of the data files

	} else {
		// Standard algorithm based on iterating the pipeline
		this->PipelineIterationAlgorithm(request,input,output);
	}


//	/* INITIALIZATION PHASE
//
//		Obtain the names of the variables that need to be time-averaged by reading
//		which variables have been loaded as cell arrays. Note that this filter cannot
//		deal with point array variables.
//
//		Recover (and broadcast) the file name of the OGS master file.
//
//		Initialize arrays for accumulating.
//
//	*/
//	// Define a vector containing the names of the arrays to compute and to exclude
//	std::string FileName;
//	std::vector<std::string> arrNames, arrNamesNotComputed{"e1","e2","e3",
//		"basins mask","coast mask","land mask","Okubo-Weiss mask","Q-criterion mask"};
//
//	VTKARRAY *vtkArray;
//	std::vector<field::Field<FLDARRAY>> arrFields;
//	
//	#ifdef PARAVIEW_USE_MPI
//	
//	// Thread 0 has all the information in input and output, therefore is the
//	// only computing the array which will later be broadcasted to all ranks
//	if (this->procId == 0) {
//		InitializeStatistics(arrNames,arrNamesNotComputed,input,output);
//		RecoverMasterFileName(FileName, input);
//		// Initalize accumulating fields
//		for (std::string arrName : arrNames) {
//			vtkArray = VTKARRAY::SafeDownCast( input->GetCellData()->GetArray(arrName.c_str()) );
//			arrFields.push_back( field::Field<FLDARRAY>(vtkArray->GetNumberOfTuples(),vtkArray->GetNumberOfComponents(),0.) );
//		}
//	}
//	
//	// Broadcast the information to all other threads if applicable
//	if (this->nProcs > 1) {
//		// Broadcast master file name
//		int str_len = FileName.length();
//		this->Controller->Broadcast(&str_len,1,0);
//		char buff[128] = ""; sprintf(buff,"%s",FileName.c_str());
//		this->Controller->Broadcast(buff,str_len,0);
//		FileName = std::string(buff);
//		
//		// Broadcast number of arrays to compute
//		int n_arrays = arrNames.size();
//		this->Controller->Broadcast(&n_arrays,1,0);
//		
//		// Allocate array on the other ranks
//		if (this->procId > 0) arrNames.resize(n_arrays,std::string(""));
//		for (int varId = 0; varId < arrNames.size(); ++varId) {
//			// Array length
//			str_len = arrNames[varId].length();
//			this->Controller->Broadcast(&str_len,1,0);
//			// Array name
//			std::memset(buff,0,128*sizeof(char));
//			sprintf(buff,"%s",arrNames[varId].c_str());
//			this->Controller->Broadcast(buff,str_len,0);
//			arrNames[varId] = std::string(buff);
//			// Array field
//			int n = 0, m = 0;
//			vtkArray = VTKARRAY::SafeDownCast( input->GetCellData()->GetArray(arrNames[varId].c_str()) );
//			if (vtkArray) {
//				n = vtkArray->GetNumberOfTuples();
//				m = vtkArray->GetNumberOfComponents();
//			}
//			this->Controller->Broadcast(&n,1,0);
//			this->Controller->Broadcast(&m,1,0);
//			if (this->procId > 0) 
//				arrFields.push_back( field::Field<FLDARRAY>(n,m,0.) );
//		}
//	}
//	
//	#else
//	
//	// This is the normal non-parallel algorithm
//	InitializeStatistics(arrNames,arrNamesNotComputed,input,output);
//	RecoverMasterFileName(FileName, input);
//
//	// Initalize accumulating fields
//	for (std::string arrName : arrNames) {
//		vtkArray = VTKARRAY::SafeDownCast( input->GetCellData()->GetArray(arrName.c_str()) );
//		arrFields.push_back( field::Field<FLDARRAY>(vtkArray->GetNumberOfTuples(),vtkArray->GetNumberOfComponents(),0.) );
//	}
//
//	#endif
//
//	// Read the OGS file to be able to generate the paths to the variable files
//	ogs::OGS ogsdata(FileName);
//	if (ogsdata.readMainFile() == -1) { 
//		vtkErrorMacro("Cannot read <"<<FileName<<">!\nAborting");
//		return 0;
//	}
//
//	this->UpdateProgress(0.1);
//
//	/* ACCUMULATE PHASE
//
//		Open and read the NetCDF files containing the variables that we previously
//		defined. Accumulate in the arrays from ii_start to ii_end.
//
//	*/
//	
//	int time_interval[2] = {this->ii_start,this->ii_end};
//
//	// Parallel partition
//	#ifdef PARAVIEW_USE_MPI
//
//	if (this->nProcs > 1) {
//		// Main thread (0) contains values for ii_start and ii_end
//		// They must be sent to the other processes
//		this->Controller->Broadcast(time_interval,2,0);
//
//		// Now everyone should have the range of start and end times
//		// We must split equally among the threads. We must also control
//		// that the number of threads is less or equal than the intervals
//		// requested.
//		int range = time_interval[1] - time_interval[0];
//		if (this->nProcs < range) {
//			// We split normally among processes assuming no remainder
//			int rangePerProcess = std::floor(range/this->nProcs);
//			this->ii_start = time_interval[0] + this->procId*rangePerProcess;
//			this->ii_end   = this->ii_start + rangePerProcess;
//			// Handle the remainder
//			int remainder = range - rangePerProcess*this->nProcs;
//			if (remainder > this->procId){
//				this->ii_start += this->procId;
//				this->ii_end   += this->procId + 1;
//			} else {
//				this->ii_start += remainder;
//				this->ii_end   += remainder;
//			}
//		} else {
//			// Each process will forcefully conduct one instant.
//			this->ii_start = (this->procId < time_interval[1]) ? this->procId     : time_interval[1];
//			this->ii_end   = (this->procId < time_interval[1]) ? this->procId + 1 : time_interval[1];
//		}
//	}
//
//	#endif
//
//	// Ensure the validity of the time range
//	if (time_interval[0] > time_interval[1]) {
//		if (this->procId == 0) 
//			vtkErrorMacro("End time is greater than initial time! Please select a valid time range.");
//		return 0;
//	}
//
//	// Loop the instants, prepare arrays and variables for the computation of the mean
//	field::Field<FLDARRAY> arrayTemp;
//	
//	// For each timestep
//	for(int ii = this->ii_start; ii < this->ii_end; ++ii) {
//		FLDARRAY ii_range = (FLDARRAY)(ii + 1. - this->ii_start);
//		// For each variable 
//		for (int varId = 0; varId < arrNames.size(); ++varId) {
//			if (arrNames[varId] == std::string("Velocity")) {
//				// Load the variable on a temporal field
//				arrayTemp.set_dim(arrFields[varId].get_n(),arrFields[varId].get_m());
//
//				std::vector<std::string> vel_vars(3);
//				vel_vars[0] = std::string("vozocrtx");
//				vel_vars[1] = std::string("vomecrty");
//				vel_vars[2] = std::string("vovecrtz");
//
//				if ( NetCDF::readNetCDF(ogsdata.var_path(arrNames[varId],ii).c_str(), vel_vars.data(), arrayTemp) != NETCDF_OK ) {					    								  
//					vtkErrorMacro("Cannot read variable <"<<arrNames[varId]<<"> in NetCDF! Aborting!"); return 0;
//				}
//				// Projecting the velocity field to the UVW grid is done a the end
//			} else {
//				// Load the variable on a temporal field
//				arrayTemp.set_dim(arrFields[varId].get_n(),arrFields[varId].get_m());
//
//				if ( NetCDF::readNetCDF(ogsdata.var_path(arrNames[varId],ii).c_str(), 
//					ogsdata.var_vname(arrNames[varId]), arrayTemp) != NETCDF_OK ) {
//					vtkErrorMacro("Cannot read variable <"<<arrNames[varId]<<"> in NetCDF! Aborting!"); return 0;
//				}
//			}
//			// Accumulate, implement 1 pass averaging algorithm
//			#pragma omp parallel shared(arrFields,arrayTemp) firstprivate(ii_range)
//			{
//			for (int nId = OMP_THREAD_NUM; nId < arrFields[varId].get_n(); nId += OMP_NUM_THREADS) {
//				for (int mId = 0; mId < arrFields[varId].get_m(); ++mId) { 
//					arrFields[varId][nId][mId] += (1./ii_range)*(arrayTemp[nId][mId] - arrFields[varId][nId][mId]);
//				}
//			}
//			}
//			arrayTemp.clear();
//		}
//		this->UpdateProgress(0.1+0.7/(this->ii_end-this->ii_start)*(ii-this->ii_start));
//	}
//
//	/* REDUCTION PHASE
//
//		If run in parallel with more than one rank, 
//		the reduction of the fields is done here.
//
//	*/
//	#ifdef PARAVIEW_USE_MPI
//
//	// Only reduce if we have more than 1 process
//	if (this->nProcs > 1) {
//		FLDARRAY n_range = (FLDARRAY)(this->ii_end-this->ii_start);
//		// For each variable 
//		for (int varId = 0; varId < arrNames.size(); ++varId) {
//			arrayTemp.set_dim(arrFields[varId].get_n(),arrFields[varId].get_m());
//
//			// Reduce to arrayTemp on rank 0
//			arrFields[varId] *= n_range;
//			this->Controller->Reduce(arrFields[varId].data(),arrayTemp.data(),arrayTemp.get_sz(),
//				vtkCommunicator::StandardOperations::SUM_OP,0);
//
//			// Average among processes
//			if (this->procId == 0) {
//				#pragma omp parallel shared(arrFields,arrayTemp) firstprivate(time_interval)
//				{
//				for (int nId = OMP_THREAD_NUM; nId < arrFields[varId].get_n(); nId += OMP_NUM_THREADS) {
//					for (int mId = 0; mId < arrayTemp.get_m(); ++mId) { 
//						arrFields[varId][nId][mId] = arrayTemp[nId][mId]/( (FLDARRAY)(time_interval[1] - time_interval[0]) );
//					}
//				}
//				}
//			}
//
//			arrayTemp.clear();
//		}
//	}
//
//	#endif
//
//	this->UpdateProgress(0.9);
//
//	/* FINALIZE
//
//		Finalization, the master process stores the arrays
//		inside the output.
//
//	*/
//
//	// Add arrays to output
//	if (this->procId == 0) {
//		for (int varId = 0; varId < arrNames.size(); ++varId) {
//			vtkArray = VTK::createVTKfromField<VTKARRAY,FLDARRAY>(arrNames[varId],arrFields[varId]);
//			output->GetCellData()->AddArray(vtkArray);
//			vtkArray->Delete();
//		}
//	}
//
	this->UpdateProgress(1.);
	return 1;
}

//----------------------------------------------------------------------------
void vtkOGSTimeAverage::PipelineIterationAlgorithm(vtkInformation *request, 
	vtkDataSet *inData, vtkDataSet *outData) {

	// Initialize.
	if (this->CurrentTimeIndex == 0) {
		outData->CopyStructure(inData);
		outData->GetFieldData()->PassData(inData->GetFieldData());
		outData->GetCellData()->PassData(inData->GetCellData());
		outData->GetPointData()->PassData(inData->GetPointData());
	}
	this->UpdateProgress(.25);

	// Accumulate new data.
	double weight = this->weights[this->CurrentTimeIndex];

	// Accumulate CellData
	for (int varId=0; varId < outData->GetCellData()->GetNumberOfArrays(); ++varId) {
		// Recover the array and the array name
		vtkDataArray *vtkDArray;
		vtkDArray = outData->GetCellData()->GetArray(varId);
		std::string arrName = vtkDArray->GetName();

		// Do not work with the basins, coasts mask, e1, e2 or e3
		if (std::string("e1")    == arrName)           continue;
		if (std::string("e2")    == arrName)           continue;
		if (std::string("e3")    == arrName)           continue;
		if (arrName.find("mask") != std::string::npos) continue; // a mask has been found

		// Recover Array values
		VTKARRAY *vtkInArray, *vtkOutArray;
		vtkInArray = VTKARRAY::SafeDownCast( inData->GetCellData()->GetArray(arrName.c_str()) );
		vtkOutArray = VTKARRAY::SafeDownCast( vtkDArray );

		field::Field<FLDARRAY> inArray, outArray;
		inArray  = VTK::createFieldfromVTK<VTKARRAY,FLDARRAY>(vtkInArray);
		outArray = VTK::createFieldfromVTK<VTKARRAY,FLDARRAY>(vtkOutArray);

		// Loop the mesh
		#pragma omp parallel
		{
		// Initialize to zero
		if (this->CurrentTimeIndex == 0) {
			for (int ii=OMP_THREAD_NUM; ii<outArray.get_n(); ii+=OMP_NUM_THREADS) {
				for (int jj=0; jj<outArray.get_m(); ++jj)
					outArray[ii][jj] = 0.;
			}
		}
		// Synchronize
		#pragma omp barrier 
		// Accumulate
		for (int ii=OMP_THREAD_NUM; ii<outArray.get_n(); ii+=OMP_NUM_THREADS)
			for (int jj=0; jj<outArray.get_m(); ++jj)
				outArray[ii][jj] += weight/this->sum_weights*(inArray[ii][jj] - outArray[ii][jj]);
		
		}

		// Set array back on output
		vtkOutArray = VTK::createVTKfromField<VTKARRAY,FLDARRAY>(arrName,outArray);
		outData->GetCellData()->RemoveArray(arrName.c_str());
		outData->GetCellData()->AddArray(vtkOutArray);
		
		this->UpdateProgress(.25 + .25*varId/outData->GetCellData()->GetNumberOfArrays());
	}
	this->UpdateProgress(0.5);

	// Accumulate PointData
	for (int varId=0; varId < outData->GetPointData()->GetNumberOfArrays(); ++varId) {
		// Recover the array and the array name
		vtkDataArray *vtkDArray;
		vtkDArray = outData->GetPointData()->GetArray(varId);
		std::string arrName = vtkDArray->GetName();

		// Do not work with the basins, coasts mask, e1, e2 or e3
		if (std::string("e1")    == arrName)           continue;
		if (std::string("e2")    == arrName)           continue;
		if (std::string("e3")    == arrName)           continue;
		if (arrName.find("mask") != std::string::npos) continue; // a mask has been found

		// Recover Array values
		VTKARRAY *vtkInArray, *vtkOutArray;
		vtkInArray = VTKARRAY::SafeDownCast( inData->GetPointData()->GetArray(arrName.c_str()) );
		vtkOutArray = VTKARRAY::SafeDownCast( vtkDArray );

		field::Field<FLDARRAY> inArray, outArray;
		inArray  = VTK::createFieldfromVTK<VTKARRAY,FLDARRAY>(vtkInArray);
		outArray = VTK::createFieldfromVTK<VTKARRAY,FLDARRAY>(vtkOutArray);

		// Loop the mesh
		#pragma omp parallel
		{
		// Initialize to zero
		if (this->CurrentTimeIndex == 0) {
			for (int ii=OMP_THREAD_NUM; ii<outArray.get_n(); ii+=OMP_NUM_THREADS) {
				for (int jj=0; jj<outArray.get_m(); ++jj)
					outArray[ii][jj] = 0.;
			}
		}
		// Synchronize
		#pragma omp barrier 
		// Accumulate
		for (int ii=OMP_THREAD_NUM; ii<outArray.get_n(); ii+=OMP_NUM_THREADS)
			for (int jj=0; jj<outArray.get_m(); ++jj)
				outArray[ii][jj] += weight/this->sum_weights*(inArray[ii][jj] - outArray[ii][jj]);
		
		}

		// Set array back on output
		vtkOutArray = VTK::createVTKfromField<VTKARRAY,FLDARRAY>(arrName,outArray);
		outData->GetPointData()->RemoveArray(arrName.c_str());
		outData->GetPointData()->AddArray(vtkOutArray);

		this->UpdateProgress(.5 + .25*varId/outData->GetCellData()->GetNumberOfArrays());
	}
	this->UpdateProgress(0.75);

	// Increment the time-step
//	printf("Time: <%d,%d> %s\n",this->CurrentTimeIndex,this->instants[this->CurrentTimeIndex],
//		this->TL[this->instants[this->CurrentTimeIndex]].as_string("%Y-%m-%d %H:%M:%S").c_str());
	this->CurrentTimeIndex++;

	// Continue executing or finish
	if (this->CurrentTimeIndex < this->instants.size()) {
		// There is still more to do.
		request->Set(vtkStreamingDemandDrivenPipeline::CONTINUE_EXECUTING(), 1);
	} else {
		// We are done. Finish up.
		request->Remove(vtkStreamingDemandDrivenPipeline::CONTINUE_EXECUTING());
		this->CurrentTimeIndex = 0;
  }	
}

//----------------------------------------------------------------------------
void vtkOGSTimeAverage::SetStartTI(const char *tstep) {

	// Set the start point of the TimeInterval

	Time::TimeObject TO(tstep,"%Y%m%d-%H:%M:%S");
	this->TI.set_start_time(TO);
	this->Modified();
}

void vtkOGSTimeAverage::SetEndTI(const char *tstep) {

	// Set the end point of the TimeInterval

	Time::TimeObject TO(tstep,"%Y%m%d-%H:%M:%S");
	this->TI.set_end_time(TO);
	this->Modified();
}

//----------------------------------------------------------------------------
vtkStringArray *vtkOGSTimeAverage::GetTimeValues() {

	// Recover time values as a string from the TimeList

	vtkStringArray *TimeValues;
	TimeValues = VTK::createVTKstrf("TimeValues",this->TL.len(),NULL);

	// Loop the TimeList to recover the instants
	for (int ii=0; ii<this->TL.len(); ++ii)
		TimeValues->SetValue(ii,this->TL[ii].as_string("%Y%m%d-%H:%M:%S"));

	return TimeValues;
}
