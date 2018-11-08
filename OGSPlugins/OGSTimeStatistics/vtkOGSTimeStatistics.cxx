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
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#include "vtkObjectFactory.h"

#include <string>
#include <vector>

vtkStandardNewMacro(vtkOGSTimeStatistics);


//----------------------------------------------------------------------------
vtkOGSTimeStatistics::vtkOGSTimeStatistics() {
	this->TimeValues = vtkStringArray::New();

	this->start_time   = 0;
	this->end_time     = 0;
	this->current_time = 0;

	this->abort = 0;
}

//----------------------------------------------------------------------------
vtkOGSTimeStatistics::~vtkOGSTimeStatistics() {
	this->TimeValues->Delete();
}

//----------------------------------------------------------------------------
int vtkOGSTimeStatistics::RequestData(vtkInformation *request, 
	vtkInformationVector **inputVector, vtkInformationVector *outputVector) {

	// get the info objects
	vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
	vtkInformation *outInfo = outputVector->GetInformationObject(0);

	// get the input and output
	vtkDataSet *input = vtkDataSet::SafeDownCast(
		inInfo->Get(vtkDataObject::DATA_OBJECT()));
	vtkDataSet *output = vtkDataSet::SafeDownCast(
		outInfo->Get(vtkDataObject::DATA_OBJECT()));

	// Check that we can perform the time statistics
	if (this->start_time > this->end_time) {
		vtkErrorMacro("Start time <"<<this->TimeValues->GetValue(this->start_time)<<
			"> is bigger than the end time <"<<this->TimeValues->GetValue(this->end_time)<<
			">.\nPlease input an end time that is bigger than the start time!");
		return 0;
	}
	if (this->start_time == this->end_time) {
		vtkErrorMacro("Start time <"<<this->TimeValues->GetValue(this->start_time)<<
			"> is equal to the end time <"<<this->TimeValues->GetValue(this->end_time)<<
			">.\nIt doesn't make sense to perform temporal statistics!");
		return 0;
	}
	if (this->abort) return 0;

	// Now we should be able to perform statistics.
	// The current time should be mapped to the start time,
	// which is smaller than the end time.

	if (this->current_time == this->start_time) 
		// This is the first iteration, we should initialize the
		// statistics.
		this->InitializeStatistics(input,output);
	else
		// We should accumulate to the statistics
		this->AccumulateStatistics(input,output);

	// Proceed to the next timestep
	this->current_time++;

	// Do we have more work to do or can we stop?
	if (this->current_time > this->end_time){
		// We are finished
		this->FinalizeStatistics(input,output);
		request->Remove(vtkStreamingDemandDrivenPipeline::CONTINUE_EXECUTING());
		this->current_time = 0;
	} else {
		// There is still more to do
		request->Set(vtkStreamingDemandDrivenPipeline::CONTINUE_EXECUTING(), 1);
	}

	// Update progress and leave
	this->UpdateProgress(1.);
	return 1;
}

//----------------------------------------------------------------------------
int vtkOGSTimeStatistics::RequestInformation(vtkInformation* vtkNotUsed(request),
  vtkInformationVector **inputVector, vtkInformationVector* outputVector) {

  	// Grab the input
  	vtkInformation *inInfo  = inputVector[0]->GetInformationObject(0);
	vtkInformation *outInfo = outputVector->GetInformationObject(0);

  	vtkDataSet *input = vtkDataSet::SafeDownCast(
		inInfo->Get(vtkDataObject::DATA_OBJECT()));

  	// Recover datevec
  	vtkStringArray *vtkdatevec = vtkStringArray::SafeDownCast(
  		input->GetFieldData()->GetAbstractArray("Datevec"));

  	if (!vtkdatevec) {
  		vtkErrorMacro("Datevec field array must be loaded!");
  		int abort = 1;
  		return 0;
  	}

  	if (this->TimeValues) {
  		this->TimeValues->Delete();
  		this->TimeValues = vtkStringArray::New();
  	}

  	// Set the time data array
  	for(int ii = 0; ii < vtkdatevec->GetNumberOfTuples(); ii++)
  		this->TimeValues->InsertNextValue(vtkdatevec->GetValue(ii));

  	// The output data of this filter has no time associated with it.  It is the
  	// result of computations that happen over all time.
  	outInfo->Remove(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
  	outInfo->Remove(vtkStreamingDemandDrivenPipeline::TIME_RANGE());

	return 1;
}

//-----------------------------------------------------------------------------
int vtkOGSTimeStatistics::RequestUpdateExtent(vtkInformation *vtkNotUsed(request),
	vtkInformationVector **inputVector, vtkInformationVector *vtkNotUsed(outputVector)) {

	vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);

	// The RequestData method will tell the pipeline executive to iterate the
	// upstream pipeline to get each time step in order.  The executive in turn
	// will call this method to get the extent request for each iteration (in this
	// case the time step).
	double *inTimes = inInfo->Get(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
	if (inTimes)
		inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP(), inTimes[this->current_time]);

	return 1;
}

//----------------------------------------------------------------------------
void vtkOGSTimeStatistics::InitializeStatistics(vtkDataSet *input, vtkDataSet *output) {
	output->CopyStructure(input);

	// Copy cell arrays
	int nCArr = input->GetCellData()->GetNumberOfArrays();
	for (int arrId = 0; arrId < nCArr; arrId++)
		output->GetCellData()->AddArray( 
			vtkFloatArray::SafeDownCast(input->GetCellData()->GetArray(arrId))
			);

	// Copy point arrays
	int nPArr = input->GetPointData()->GetNumberOfArrays();
	for (int arrId = 0; arrId < nPArr; arrId++)
		output->GetPointData()->AddArray(
			vtkFloatArray::SafeDownCast(input->GetPointData()->GetArray(arrId))
			);

	// Copy field arrays
	int nFArr = input->GetFieldData()->GetNumberOfArrays();
	for (int arrId = 0; arrId < nFArr; arrId++)
		output->GetFieldData()->AddArray(
			vtkStringArray::SafeDownCast(input->GetFieldData()->GetAbstractArray(arrId))
			);
}

//----------------------------------------------------------------------------
void vtkOGSTimeStatistics::AccumulateStatistics(vtkDataSet *input, vtkDataSet *output) {

	// Recover number of cells and cell arrays
	int nCells   = input->GetNumberOfCells();
	int nCellArr = input->GetCellData()->GetNumberOfArrays();
	// Recover number of points and point arrays
	int nP    = input->GetNumberOfPoints();
	int nPArr = input->GetPointData()->GetNumberOfArrays();

	// Recover cell and point arrays
	std::vector<vtkFloatArray*> vtkICellArrays,  vtkOCellArrays;
	std::vector<vtkFloatArray*> vtkIPointArrays, vtkOPointArrays;
	for(int arrId = 0; arrId < nCellArr; arrId++) {		
		vtkICellArrays.push_back(
			vtkFloatArray::SafeDownCast(input->GetCellData()->GetArray(arrId))
				);
		vtkOCellArrays.push_back(
			vtkFloatArray::SafeDownCast(output->GetCellData()->GetArray(arrId))
				);
	}
	for(int arrId = 0; arrId < nPArr; arrId++) {
		vtkIPointArrays.push_back(
			vtkFloatArray::SafeDownCast(input->GetPointData()->GetArray(arrId))
				);
		vtkOPointArrays.push_back(
			vtkFloatArray::SafeDownCast(output->GetPointData()->GetArray(arrId))
				);
	}

	// Loop the mesh, to save iterations get the maximum of the cells and points
	for (int id = 0; id < std::max(nCells,nP); id++) {
		// Accumulate Cell Arrays
		for(int arrId = 0; id < nCells && arrId < nCellArr; arrId++) {
			int ncomp = vtkICellArrays[arrId]->GetNumberOfComponents();
			// Scalar Array
			if (ncomp == 1) {
				double valI = vtkICellArrays[arrId]->GetTuple1(id);
				double valO = vtkOCellArrays[arrId]->GetTuple1(id);
				// Set Value
				vtkOCellArrays[arrId]->SetTuple1(id,valI+valO);
			}
			// Vector array
			if (ncomp == 3) {
				double valI[3]; vtkICellArrays[arrId]->GetTuple(id,valI);
				double valO[3]; vtkOCellArrays[arrId]->GetTuple(id,valO);
				// Set Value
				for (int ii=0; ii<3; ii++) valO[ii] += valI[ii];
				vtkOCellArrays[arrId]->SetTuple(id,valO);
			}
		}
		// Accumulate Cell Arrays
		for(int arrId = 0; id < nP && arrId < nPArr; arrId++) {
			int ncomp = vtkIPointArrays[arrId]->GetNumberOfComponents();
			// Scalar Array
			if (ncomp == 1) {
				double valI = vtkIPointArrays[arrId]->GetTuple1(id);
				double valO = vtkOPointArrays[arrId]->GetTuple1(id);
				// Set Value
				vtkOPointArrays[arrId]->SetTuple1(id,valI+valO);
			}
			// Vector array
			if (ncomp == 3) {
				double valI[3]; vtkIPointArrays[arrId]->GetTuple(id,valI);
				double valO[3]; vtkOPointArrays[arrId]->GetTuple(id,valO);
				// Set Value
				for (int ii=0; ii<3; ii++) valO[ii] += valI[ii];
				vtkOPointArrays[arrId]->SetTuple(id,valO);
			}
		}
		// Update progress
		this->UpdateProgress(0.+0.9/(std::max(nCells,nP))*id);
	}
}

//----------------------------------------------------------------------------
void vtkOGSTimeStatistics::FinalizeStatistics(vtkDataSet *input, vtkDataSet *output) {

	// Recover number of cells and cell arrays
	int nCells   = input->GetNumberOfCells();
	int nCellArr = input->GetCellData()->GetNumberOfArrays();
	// Recover number of points and point arrays
	int nP    = input->GetNumberOfPoints();
	int nPArr = input->GetPointData()->GetNumberOfArrays();

	// Current time is also used as a counter
	double divi = (double)(this->current_time);

	// Recover cell and point arrays
	std::vector<vtkFloatArray*> vtkOCellArrays,vtkOPointArrays;
	for(int arrId = 0; arrId < nCellArr; arrId++) {		
		vtkOCellArrays.push_back(
			vtkFloatArray::SafeDownCast(output->GetCellData()->GetArray(arrId))
				);
	}
	for(int arrId = 0; arrId < nPArr; arrId++) {
		vtkOPointArrays.push_back(
			vtkFloatArray::SafeDownCast(output->GetPointData()->GetArray(arrId))
				);
	}

	// Loop the mesh, to save iterations get the maximum of the cells and points
	for (int id = 0; id < std::max(nCells,nP); id++) {
		// Accumulate Cell Arrays
		for(int arrId = 0; id < nCells && arrId < nCellArr; arrId++) {
			int ncomp = vtkOCellArrays[arrId]->GetNumberOfComponents();
			// Scalar Array
			if (ncomp == 1) {
				double valO = vtkOCellArrays[arrId]->GetTuple1(id);
				// Set Value
				vtkOCellArrays[arrId]->SetTuple1(id,valO/divi);
			}
			// Vector array
			if (ncomp == 3) {
				double valO[3]; vtkOCellArrays[arrId]->GetTuple(id,valO);
				// Set Value
				for (int ii=0; ii<3; ii++) valO[ii] /= divi;
				vtkOCellArrays[arrId]->SetTuple(id,valO);
			}
		}
		// Accumulate Cell Arrays
		for(int arrId = 0; id < nP && arrId < nPArr; arrId++) {
			int ncomp = vtkOPointArrays[arrId]->GetNumberOfComponents();
			// Scalar Array
			if (ncomp == 1) {
				double valO = vtkOPointArrays[arrId]->GetTuple1(id);
				// Set Value
				vtkOPointArrays[arrId]->SetTuple1(id,valO/divi);
			}
			// Vector array
			if (ncomp == 3) {
				double valO[3]; vtkOPointArrays[arrId]->GetTuple(id,valO);
				// Set Value
				for (int ii=0; ii<3; ii++) valO[ii] /= divi;
				vtkOPointArrays[arrId]->SetTuple(id,valO);
			}
		}
		// Update progress
		this->UpdateProgress(0.9+0.1/(std::max(nCells,nP))*id);
	}
}

//----------------------------------------------------------------------------
void vtkOGSTimeStatistics::SetStartTime(const char *tstep) {
	// Obtain the timestep index
	for (int ii = 0; ii < this->TimeValues->GetNumberOfTuples(); ii++)
		if (std::string(this->TimeValues->GetValue(ii)) == std::string(tstep))
			this->start_time = ii;
	// Rewind to start
	this->current_time = this->start_time;
	this->Modified();
}

void vtkOGSTimeStatistics::SetEndTime(const char *tstep) {
	// Obtain the timestep index
	for (int ii = 0; ii < this->TimeValues->GetNumberOfTuples(); ii++)
		if (std::string(this->TimeValues->GetValue(ii)) == std::string(tstep))
			this->end_time = ii;
	this->Modified();
}

//----------------------------------------------------------------------------
vtkStringArray *vtkOGSTimeStatistics::GetTimeValues() {
	return this->TimeValues;
}
