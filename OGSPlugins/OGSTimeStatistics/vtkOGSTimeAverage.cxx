/*=========================================================================

  Program:   OGSTimeAverage
  Module:    vtkOGSTimeAverage.cxx

  Copyright (c) 2020 Arnau Miro, OGS
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
#include "netcdfio.h"
#include "OGS.hpp"
#include "fieldOperations.h"
#include "vtkFields.h"
#include "vtkOGSTimeCommons.hpp"

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
			this->instants.clear(); this->weights.clear(); // make sure vectors are empty

			if (this->TL.select(&req,this->instants,this->weights) == TIME_ERR) {
				vtkErrorMacro("Error selecting instants with Requestor! Cannot continue...");
				return 0;
			}
//			for (int ii=0; ii<instants.size();++ii)
//				printf("ind: %d, weight: %f, %s\n",this->instants[ii],this->weights[ii],
//					this->TL[instants[ii]].as_string("%Y-%m-%d %H:%M:%S").c_str());
		}

		// Broadcast instants and weights to all the processors if needed
		#ifdef PARAVIEW_USE_MPI
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
		#endif

		// Initialize sum of weights for the average
		this->sum_weights = 0.;
	}

	int retval = 0;
	if (this->use_files) { 
		// Optimized algorithm based on parallel access of the data files
		retval = this->FileIterationAlgorithm(request,input,output);
	} else {
		// Standard algorithm based on iterating the pipeline
		retval = this->PipelineIterationAlgorithm(request,input,output);
	}

	this->UpdateProgress(1.);
	return retval;
}

//----------------------------------------------------------------------------
int vtkOGSTimeAverage::PipelineIterationAlgorithm(vtkInformation *request, 
	vtkDataSet *inData, vtkDataSet *outData) {

	// Stop all threads except from the master to execute
	#ifdef PARAVIEW_USE_MPI
	if (this->procId > 0) return 1;
	#endif

	// Initialize.
	if (this->CurrentTimeIndex == 0) {
		outData->CopyStructure(inData);
		outData->GetFieldData()->PassData(inData->GetFieldData());
		outData->GetCellData()->PassData(inData->GetCellData());
		outData->GetPointData()->PassData(inData->GetPointData());
	}
	this->UpdateProgress(.25);

	// Accumulate new data.
	int    itime  = this->instants[this->CurrentTimeIndex];
	double weight = this->weights[this->CurrentTimeIndex];
	this->sum_weights += weight;
//	printf("Time: <%d,%d> %s\n",this->CurrentTimeIndex,itime,this->TL[itime].as_string("%Y-%m-%d %H:%M:%S").c_str());

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
		if (this->CurrentTimeIndex == 0) outArray.set_val(0.);

		// Loop the mesh
		#pragma omp parallel shared(inArray,outArray)
		{
		for (int ii=OMP_THREAD_NUM; ii<outArray.get_n(); ii+=OMP_NUM_THREADS)
			for (int jj=0; jj<outArray.get_m(); ++jj)
				outArray[ii][jj] += weight/this->sum_weights*(inArray[ii][jj] - outArray[ii][jj]);
		}

		// Set array back on output
		vtkOutArray = VTK::createVTKfromField<VTKARRAY,FLDARRAY>(arrName,outArray);
		outData->GetCellData()->RemoveArray(arrName.c_str());
		outData->GetCellData()->AddArray(vtkOutArray);
		vtkOutArray->Delete();
		
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
		if (this->CurrentTimeIndex == 0) outArray.set_val(0.);

		// Loop the mesh
		#pragma omp parallel shared(inArray,outArray)
		{
		// Accumulate
		for (int ii=OMP_THREAD_NUM; ii<outArray.get_n(); ii+=OMP_NUM_THREADS)
			for (int jj=0; jj<outArray.get_m(); ++jj)
				outArray[ii][jj] += weight/this->sum_weights*(inArray[ii][jj] - outArray[ii][jj]);
		}

		// Set array back on output
		vtkOutArray = VTK::createVTKfromField<VTKARRAY,FLDARRAY>(arrName,outArray);
		outData->GetPointData()->RemoveArray(arrName.c_str());
		outData->GetPointData()->AddArray(vtkOutArray);
		vtkOutArray->Delete();

		this->UpdateProgress(.5 + .25*varId/outData->GetPointData()->GetNumberOfArrays());
	}
	this->UpdateProgress(0.75);

	// Increment the time-step
	this->CurrentTimeIndex++;

	// Continue executing or finish
	if (this->CurrentTimeIndex < this->instants.size()) {
		// There is still more to do.
		request->Set(vtkStreamingDemandDrivenPipeline::CONTINUE_EXECUTING(), 1);
	} else {
		// We are done. Finish up.
		request->Remove(vtkStreamingDemandDrivenPipeline::CONTINUE_EXECUTING());
		this->CurrentTimeIndex = 0;
		this->sum_weights = 0.;
	}

	return 1;
}

//----------------------------------------------------------------------------
int vtkOGSTimeAverage::FileIterationAlgorithm(vtkInformation *request, 
	vtkDataSet *inData, vtkDataSet *outData) {

	// Check inputs
	if (!inData->IsA("vtkRectilinearGrid")) {
		vtkErrorMacro("A rectilinear grid is needed for this method! Cannot continue...");
		return 0;		
	}

	if (inData->GetPointData()->GetNumberOfArrays() > 0) {
		vtkErrorMacro("Data contains point arrays! Cannot continue...");
		return 0;
	}

	/* INITIALIZATION PHASE

		Obtain the names of the variables that need to be time-averaged by reading
		which variables have been loaded as cell arrays. Note that this filter cannot
		deal with point array variables.

		Recover (and broadcast) the file name of the OGS master file.

		Initialize arrays for accumulating.

	*/
	// Define a vector containing the names of the arrays to compute and to exclude
	VTKARRAY *vtkOutArray;
	std::string FileName;
	
	#ifdef PARAVIEW_USE_MPI
	
	// Thread 0 has all the information in input and output, therefore is the
	// only computing the array which will later be broadcasted to all ranks
	if (this->procId == 0) {
		// Initialize statistics
		outData->CopyStructure(inData);
		outData->GetFieldData()->PassData(inData->GetFieldData());
		outData->GetCellData()->PassData(inData->GetCellData());

		RecoverMasterFileName(FileName, inData);
	}
	
	// Broadcast the information to all other threads if applicable
	if (this->nProcs > 1) {
		// Broadcast master file name
		char buff[128] = ""; sprintf(buff,"%s",FileName.c_str());
		this->Controller->Broadcast(buff,FileName.length(),0);
		if (procId > 0) FileName = std::string(buff);

		// Broadcast output data and its arrays
		this->Controller->Broadcast(outData,0);
	}
	
	#else
	
	// This is the normal non-parallel algorithm
	outData->CopyStructure(inData);
	outData->GetFieldData()->PassData(inData->GetFieldData());
	outData->GetCellData()->PassData(inData->GetCellData());

	RecoverMasterFileName(FileName, inData);

	#endif

	// Initialize arrays
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
		
		vtkOutArray = VTKARRAY::SafeDownCast(vtkDArray);

		field::Field<FLDARRAY> outArray;
		outArray = VTK::createFieldfromVTK<VTKARRAY,FLDARRAY>(vtkOutArray);
		outArray.set_val(0.);

		// Set array back on output
		vtkOutArray = VTK::createVTKfromField<VTKARRAY,FLDARRAY>(arrName,outArray);
		outData->GetCellData()->RemoveArray(arrName.c_str());
		outData->GetCellData()->AddArray(vtkOutArray);
		vtkOutArray->Delete();
	}

	this->UpdateProgress(0.1);

	// Read the OGS file to be able to generate the paths to the variable files
	ogs::OGS ogsdata(FileName);
	if (ogsdata.readMainFile() == -1) { 
		vtkErrorMacro("Cannot read <"<<FileName<<">!\nAborting");
		return 0;
	}

	// Read mesh. ProjID is not important since only the dimensions
	// and the e1, e2 and e3 arrays are necessary.
	if (ogsdata.readMesh(0) < 0) {
		vtkErrorMacro("Problems reading the mesh!\nAborting.");
		return 0;	
	}
	
	/* ACCUMULATE PHASE

		Instants and weights have already been broadcasted
		so at this point all the processors should be able
		to move forward to the accumulation phase.

		Open and read the NetCDF files containing the variables defined
		in the outData structure.
	*/

	// Parallel partition, define starting and ending ranges for each worker
	int ii_start = 0, ii_end = this->instants.size();

	#ifdef PARAVIEW_USE_MPI
	if (this->nProcs > 1) {
		// Everyone should have the range of start and end times
		// We must split equally among the threads. We must also control
		// that the number of threads is less or equal than the intervals
		// requested.
		int range = ii_end - ii_start;
		if (this->nProcs < range) {
			// We split normally among processes assuming no remainder
			int rangePerProcess = std::floor(range/this->nProcs);
			ii_start = ii_start + this->procId*rangePerProcess;
			ii_end   = ii_start + rangePerProcess;
			// Handle the remainder
			int remainder = range - rangePerProcess*this->nProcs;
			if (remainder > this->procId){
				ii_start += this->procId;
				ii_end   += this->procId + 1;
			} else {
				ii_start += remainder;
				ii_end   += remainder;
			}
		} else {
			// Each process will forcefully conduct one instant.
			ii_start = (this->procId < ii_end) ? this->procId     : ii_end;
			ii_end   = (this->procId < ii_end) ? this->procId + 1 : ii_end;
		}
	}
	#endif

	// Loop the instants
	double w_sum = 0.;
	for(int inst = ii_start; inst < ii_end; ++inst) {
		int itime     = this->instants[inst];
		double weight = this->weights[inst];
		w_sum += weight;
//		printf("Time: <%d,%d> %s\n",inst,itime,this->TL[itime].as_string("%Y-%m-%d %H:%M:%S").c_str());
		// Loop cell arrays
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

			vtkOutArray = VTKARRAY::SafeDownCast(vtkDArray);

			field::Field<FLDARRAY> outArray, inArray;
			outArray = VTK::createFieldfromVTK<VTKARRAY,FLDARRAY>(vtkOutArray);
			inArray.set_dim(outArray.get_n(),outArray.get_m());

			// Deal with the velocity
			if (arrName == std::string("Velocity")) {

				field::Field<FLDARRAY> tmp(outArray.get_n(),outArray.get_m(),0.);

				std::vector<std::string> vel_vars;
				strsplit(ogsdata.var_vname(arrName),",",vel_vars);

				if ( NetCDF::readNetCDF(ogsdata.var_path(arrName,itime).c_str(), vel_vars.data(), tmp) != NETCDF_OK ) {
					vtkErrorMacro("Cannot read variable <"<<arrName<<"> in NetCDF! Aborting!"); 
					return 0;
				}

				// Projecting the velocity field to the UVW grid
				inArray = field::UVW2T(tmp,ogsdata.e1(),ogsdata.e2(),ogsdata.e3(),ogsdata.nlon()-1,ogsdata.nlat()-1,ogsdata.nlev()-1);
			} else {
				if ( NetCDF::readNetCDF(ogsdata.var_path(arrName,itime).c_str(),ogsdata.var_vname(arrName), inArray) != NETCDF_OK ) {
					vtkErrorMacro("Cannot read variable <"<<arrName<<" ("<<ogsdata.var_vname(arrName)<<")> in NetCDF! Aborting!"); 
					return 0;
				}
			}
			// Accumulate, implement 1 pass averaging algorithm
			#pragma omp parallel shared(inArray,outArray)
			{
			for (int ii = OMP_THREAD_NUM; ii < outArray.get_n(); ii += OMP_NUM_THREADS)
				for (int jj = 0; jj < outArray.get_m(); ++jj)
					outArray[ii][jj] += weight*inArray[ii][jj];
//					outArray[ii][jj] += (weight/w_sum)*(inArray[ii][jj] - outArray[ii][jj]);
			}
			// Set array back on output
			vtkOutArray = VTK::createVTKfromField<VTKARRAY,FLDARRAY>(arrName,outArray);
			outData->GetCellData()->RemoveArray(arrName.c_str());
			outData->GetCellData()->AddArray(vtkOutArray);
			vtkOutArray->Delete();

		}
		this->UpdateProgress(0.1+0.7/(ii_end-ii_start)*(inst-ii_start));
	}

	// At this point, outData contains the partial sums of the different arrays
	// Now we need to loop the arrays in outData and reduce them to proc 0
	// This operation needs to be done only for prodId > 0
	#ifdef PARAVIEW_USE_MPI
	if (this->nProcs > 1) {
		// Loop cell arrays
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

			this->Controller->AllReduce(vtkDArray,vtkDArray,vtkCommunicator::StandardOperations::SUM_OP);
		
			outData->GetCellData()->RemoveArray(arrName.c_str());
			outData->GetCellData()->AddArray(vtkDArray);
		}
	}
	#endif
	this->UpdateProgress(0.9);

	// Now proc 0 has all the information. We need to loop the mesh once more to finish
	// the average.
	if (this->procId == 0) {
		// Loop cell arrays
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
			
			vtkOutArray = VTKARRAY::SafeDownCast(vtkDArray);

			field::Field<FLDARRAY> outArray;
			outArray = VTK::createFieldfromVTK<VTKARRAY,FLDARRAY>(vtkOutArray);

			// Accumulate, implement 1 pass averaging algorithm
			#pragma omp parallel shared(outArray)
			{
			for (int ii = OMP_THREAD_NUM; ii < outArray.get_n(); ii += OMP_NUM_THREADS)
				for (int jj = 0; jj < outArray.get_m(); ++jj)
					outArray[ii][jj] /= w_sum;
			}

			// Set array back on output
			vtkOutArray = VTK::createVTKfromField<VTKARRAY,FLDARRAY>(arrName,outArray);
			outData->GetCellData()->RemoveArray(arrName.c_str());
			outData->GetCellData()->AddArray(vtkOutArray);
			vtkOutArray->Delete();
		}
	}
	this->UpdateProgress(0.95);

	// And we're done!
	return 1;
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
