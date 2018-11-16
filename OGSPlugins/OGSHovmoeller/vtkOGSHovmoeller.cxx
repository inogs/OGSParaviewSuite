/*=========================================================================

  Program:   OGSHovmoeller
  Module:    vtkOGSHovmoeller.cxx

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkOGSHovmoeller.h"

#include "vtkDataObject.h"
#include "vtkDataObjectTypes.h"
#include "vtkCell.h"
#include "vtkCellData.h"
#include "vtkCellLocator.h"
#include "vtkCharArray.h"
#include "vtkGenericCell.h"
#include "vtkPointData.h"
#include "vtkRectilinearGrid.h"
#include "vtkFloatArray.h"
#include "vtkStringArray.h"
#include "vtkPointSet.h"
#include "vtkDataSet.h"
#include "vtkTable.h"
#include "vtkExecutive.h"
#include "vtkSmartPointer.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkStaticCellLocator.h"
#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"

namespace VTK
{
// Include the VTK Operations
#include "../_utils/vtkOperations.cpp"
}

#define CELL_TOLERANCE_FACTOR_SQR 1e-6

vtkStandardNewMacro(vtkOGSHovmoeller);
vtkCxxSetObjectMacro(vtkOGSHovmoeller, CellLocatorPrototype, vtkAbstractCellLocator);

//----------------------------------------------------------------------------
vtkOGSHovmoeller::vtkOGSHovmoeller() {
	this->SetNumberOfInputPorts(2);
	this->SetNumberOfOutputPorts(1);

	this->PointList  = nullptr;
	this->CellList   = nullptr;
	this->field      = nullptr;

	this->CellLocatorPrototype = nullptr;

	this->TimeValues = vtkStringArray::New();

	this->start_time   = 0;
	this->end_time     = 0;
	this->current_time = 0;
	this->DepthScale   = 1000.;

	this->abort = 0;
}

//----------------------------------------------------------------------------
vtkOGSHovmoeller::~vtkOGSHovmoeller() {
	delete this->PointList;
	delete this->CellList;

	this->vtkOGSHovmoeller::SetCellLocatorPrototype(nullptr);
}

//----------------------------------------------------------------------------
int vtkOGSHovmoeller::RequestData(vtkInformation *request, 
	vtkInformationVector **inputVector, vtkInformationVector *outputVector) {

	// get the info objects
	vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
	vtkInformation *sourceInfo = inputVector[1]->GetInformationObject(0);
	vtkInformation *outInfo = outputVector->GetInformationObject(0);

	// input contains the interpolating line information (number of points, etc)
	vtkDataSet *input = vtkDataSet::SafeDownCast(
		inInfo->Get(vtkDataObject::DATA_OBJECT()));
	// Source contains the data where to interpolate from
	vtkDataSet *source = vtkDataSet::SafeDownCast(
		sourceInfo->Get(vtkDataObject::DATA_OBJECT()));
	// Output is a vtkTable with the interpolated data per each timestep
	vtkTable *output = vtkTable::SafeDownCast(
		outInfo->Get(vtkDataObject::DATA_OBJECT()));

	// Check that we can perform the Hovmoeller plot
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

	// Write the coordinates as a column on the table
	// This action is only performed at the beginning
	if (this->current_time == this->start_time) {
		int npoints = input->GetNumberOfPoints();
		vtkFloatArray *vtkDepth = VTK::createVTKscaf("depth",npoints,NULL);
		for (int pId = 0; pId < npoints; pId++) {
			double xyz[3]; input->GetPoint(pId,xyz);
			vtkDepth->SetTuple1(pId,xyz[2]/this->DepthScale);
		}
		output->AddColumn(vtkDepth);
	}

	// Now we should be able to perform statistics.
	// The current time should be mapped to the start time,
	// which is smaller than the end time.

	this->Interpolate(input,source,output);

	// Proceed to the next timestep
	this->current_time++;

	// Do we have more work to do or can we stop?
	if (this->current_time > this->end_time){
		// We are finished
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
int vtkOGSHovmoeller::RequestInformation( vtkInformation *vtkNotUsed(request),
	vtkInformationVector **inputVector, vtkInformationVector *outputVector) {

	// get the info objects
	vtkInformation *inInfo     = inputVector[0]->GetInformationObject(0);
	vtkInformation *sourceInfo = inputVector[1]->GetInformationObject(0);
	vtkInformation *outInfo    = outputVector->GetInformationObject(0);

	vtkDataSet *source = vtkDataSet::SafeDownCast(
		sourceInfo->Get(vtkDataObject::DATA_OBJECT()));

	// Recover datevec
	vtkStringArray *vtkdatevec = vtkStringArray::SafeDownCast(
		source->GetFieldData()->GetAbstractArray("Datevec"));

	if (vtkdatevec) {
		if (this->TimeValues) {
			this->TimeValues->Delete(); 		
			this->TimeValues = vtkStringArray::New();
		}
		// Set the time data array
		for(int ii = 0; ii < vtkdatevec->GetNumberOfTuples(); ii++)
			this->TimeValues->InsertNextValue(vtkdatevec->GetValue(ii));
	}

	// The output data of this filter has no time associated with it.  It is the
	// result of computations that happen over all time.
	outInfo->Remove(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
	outInfo->Remove(vtkStreamingDemandDrivenPipeline::TIME_RANGE());

	return 1;
}

//-----------------------------------------------------------------------------
int vtkOGSHovmoeller::RequestUpdateExtent(vtkInformation *vtkNotUsed(request),
	vtkInformationVector **inputVector, vtkInformationVector *vtkNotUsed(outputVector)) {

	vtkInformation *inInfo     = inputVector[0]->GetInformationObject(0);
	vtkInformation *sourceInfo = inputVector[1]->GetInformationObject(0);

	// The RequestData method will tell the pipeline executive to iterate the
	// upstream pipeline to get each time step in order.  The executive in turn
	// will call this method to get the extent request for each iteration (in this
	// case the time step).
	double *inTimes = sourceInfo->Get(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
	if (inTimes)
		sourceInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP(), inTimes[this->current_time]);

	return 1;
}

//----------------------------------------------------------------------------
int vtkOGSHovmoeller::RequestDataObject(vtkInformation* vtkNotUsed(request),
  vtkInformationVector** inputVector, vtkInformationVector* outputVector) {

	const char* outTypeStr = vtkDataObjectTypes::GetClassNameFromTypeId(VTK_TABLE);

	vtkInformation *outInfo = outputVector->GetInformationObject(0);
	vtkDataObject* output   = outInfo->Get(vtkDataObject::DATA_OBJECT());

	if (!output || !output->IsA(outTypeStr)) {
		vtkDataObject* newOutput = vtkDataObjectTypes::NewDataObject(VTK_TABLE);
		if (!newOutput) {
			vtkErrorMacro("Could not create chosen output data type: " << outTypeStr);
			return 0;
		}
		outInfo->Set(vtkDataObject::DATA_OBJECT(), newOutput);
		this->GetOutputPortInformation(0)->Set(
			vtkDataObject::DATA_EXTENT_TYPE(), newOutput->GetExtentType());
		newOutput->Delete();
	}

	return 1;
}

//----------------------------------------------------------------------------
int vtkOGSHovmoeller::FillOutputPortInformation(
	int vtkNotUsed(port), vtkInformation* info) {

	info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataObject");
	return 1;
}

//----------------------------------------------------------------------------
void vtkOGSHovmoeller::Interpolate(vtkDataSet *input, vtkDataSet *source, vtkTable *output) {

	// Build Point List
	delete this->PointList;
	this->PointList = new vtkDataSetAttributes::FieldList(1);
	this->PointList->InitializeFieldList(source->GetPointData());

	// Build Cell List
	delete this->CellList;
	this->CellList = new vtkDataSetAttributes::FieldList(1);
	this->CellList->InitializeFieldList(source->GetCellData());

	// Preallocate weights
	double *weights;
	weights = new double[source->GetMaxCellSize()];

	// Create the cell locator object
	vtkCellLocator *cellLocator = vtkCellLocator::New();
	cellLocator->SetDataSet(source); cellLocator->BuildLocator();

	// Recover the number of points in the input
	int npoints = input->GetNumberOfPoints();

	// Create auxiliary point data
	vtkPointData *auxPD = vtkPointData::New();

	// Allocate storage for output PointData
	// All input PD is passed to output as PD. Those arrays in input CD that are
	// not present in output PD will be passed as output PD.
	auxPD->InterpolateAllocate((*this->PointList), npoints, npoints);
	// We assume we either have points or cells
	auxPD->InterpolateAllocate((*this->CellList), npoints, npoints);

	// Initialize Output Arrays
	int nArr = auxPD->GetNumberOfArrays();
	for (int arrId = 0; arrId < nArr; arrId++) {
		auxPD->GetArray(arrId)->SetNumberOfTuples(npoints);
		auxPD->GetArray(arrId)->Fill(0);
	}

	// Loop the number of points
	vtkNew<vtkGenericCell> gcell;
	int abort = 0;
	for (int ii = 0; ii < npoints && !abort; ii++) {
		// Update progress
		this->UpdateProgress((double)(ii)/(double)(npoints));
		abort = GetAbortExecute();

		// Get the xyz coordinate of the point in the input dataset
		double xyz[3]; input->GetPoint(ii, xyz);

		// Find the cell that contains xyz and get it
		vtkIdType cellId; int subId;
		double dist2, pcoords[3], closestPoint[3];
		cellLocator->FindClosestPoint(xyz, closestPoint, gcell.GetPointer(),cellId,subId,dist2);

		// Evaluate interpolation weights
		vtkCell* cell = nullptr;
		if (cellId >= 0) {
			cell = source->GetCell(cellId);
			// Compute a tolerance proportional to the cell length.
			cell->EvaluatePosition(xyz, closestPoint, subId, pcoords, dist2, weights);
			// Abort if the distance is too big
			if (dist2 > (cell->GetLength2() * CELL_TOLERANCE_FACTOR_SQR))
				continue;
		}

		// Check if the cell has been found
		if (cell) {
			// Interpolate point array data
			auxPD->InterpolatePoint(
					(*this->PointList), source->GetPointData(), 0, ii, cell->PointIds, weights
				);
			// Interpolate cell array data
			for (int arrId=0; arrId < source->GetCellData()->GetNumberOfArrays(); arrId++) {
				auxPD->CopyTuple(
					source->GetCellData()->GetArray(arrId),
					auxPD->GetArray(arrId),
					cellId,ii);
			}
		}
	}

	// We can now write the variable in the table
  	vtkStringArray *vtkdate = vtkStringArray::SafeDownCast(
  		source->GetFieldData()->GetAbstractArray("Date"));

	vtkFloatArray *vtkVar =vtkFloatArray::SafeDownCast(
		auxPD->GetArray(this->field));
	vtkVar->SetName(vtkdate->GetValue(0).c_str());

	output->AddColumn(vtkVar);

	// Delete and deallocate
	delete [] weights;
	cellLocator->Delete();
	auxPD->Delete();
}

//----------------------------------------------------------------------------
void vtkOGSHovmoeller::SetSourceConnection(vtkAlgorithmOutput* algOutput) {
	this->SetInputConnection(1, algOutput);
}

void vtkOGSHovmoeller::SetSourceData(vtkDataObject *input) {
	this->SetInputData(1, input);
}

vtkDataObject *vtkOGSHovmoeller::GetSource() {
	if (this->GetNumberOfInputConnections(1) < 1)
		return nullptr;

	return this->GetExecutive()->GetInputData(1, 0);
}

//----------------------------------------------------------------------------
void vtkOGSHovmoeller::SetStartTime(const char *tstep) {
	// Obtain the timestep index
	for (int ii = 0; ii < this->TimeValues->GetNumberOfTuples(); ii++)
		if (std::string(this->TimeValues->GetValue(ii)) == std::string(tstep))
			this->start_time = ii;
	// Rewind to start
	this->current_time = this->start_time;
	this->Modified();
}

void vtkOGSHovmoeller::SetEndTime(const char *tstep) {
	// Obtain the timestep index
	for (int ii = 0; ii < this->TimeValues->GetNumberOfTuples(); ii++)
		if (std::string(this->TimeValues->GetValue(ii)) == std::string(tstep))
			this->end_time = ii;
	this->Modified();
}

vtkStringArray *vtkOGSHovmoeller::GetTimeValues() {
	return this->TimeValues;
}