/*=========================================================================

  Program:   OGSDepthProfile
  Module:    vtkOGSDepthProfile.cxx

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkOGSDepthProfile.h"

#include "vtkCell.h"
#include "vtkCellData.h"
#include "vtkCellLocator.h"
#include "vtkStringArray.h"
#include "vtkGenericCell.h"
#include "vtkPointData.h"
#include "vtkPointSet.h"
#include "vtkDataSet.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkStaticCellLocator.h"
#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#define CELL_TOLERANCE_FACTOR_SQR 1e-6

vtkStandardNewMacro(vtkOGSDepthProfile);
vtkCxxSetObjectMacro(vtkOGSDepthProfile, CellLocatorPrototype, vtkAbstractCellLocator);

//----------------------------------------------------------------------------
vtkOGSDepthProfile::vtkOGSDepthProfile() {
	this->SetNumberOfInputPorts(2);
	this->SetNumberOfOutputPorts(1);

	this->PointList  = nullptr;
	this->CellList   = nullptr;

	this->CellLocatorPrototype = nullptr;
}

//----------------------------------------------------------------------------
vtkOGSDepthProfile::~vtkOGSDepthProfile() {
	delete this->PointList;
	delete this->CellList;

	this->vtkOGSDepthProfile::SetCellLocatorPrototype(nullptr);
}

//----------------------------------------------------------------------------
int vtkOGSDepthProfile::RequestData(vtkInformation *vtkNotUsed(request), 
	vtkInformationVector **inputVector, vtkInformationVector *outputVector) {

	// get the info objects
	vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
	vtkInformation *sourceInfo = inputVector[1]->GetInformationObject(0);
	vtkInformation *outInfo = outputVector->GetInformationObject(0);

	// get the input and output
	// input contains the interpolating line information (number of points, etc)
	vtkDataSet *input = vtkDataSet::SafeDownCast(
		inInfo->Get(vtkDataObject::DATA_OBJECT()));
	// Source contains the data where to interpolate from
	vtkDataSet *source = vtkDataSet::SafeDownCast(
		sourceInfo->Get(vtkDataObject::DATA_OBJECT()));
	vtkDataSet *output = vtkDataSet::SafeDownCast(
		outInfo->Get(vtkDataObject::DATA_OBJECT()));

	// First, copy the input to the output as a starting point
	output->CopyStructure(input);
	output->GetPointData()->SetCopyAttribute(vtkDataSetAttributes::SCALARS,
		2,vtkDataSetAttributes::INTERPOLATE);

	// If there is data to interpolate, begin the interpolation
	if (source) {
		this->Initialize(input, source, output);
		this->Interpolate(input, source, output);
		// Remove some arrays
		output->GetPointData()->RemoveArray("basins mask");
		output->GetPointData()->RemoveArray("coast mask");
		output->GetPointData()->RemoveArray("e1");
		output->GetPointData()->RemoveArray("e2");
		output->GetPointData()->RemoveArray("e3");
		// Make sure the metadata array is passed to the output
		vtkStringArray *vtkmetadata = vtkStringArray::SafeDownCast(
			source->GetFieldData()->GetAbstractArray("Metadata"));
		output->GetFieldData()->AddArray(vtkmetadata);
	}

	return 1;
}

//----------------------------------------------------------------------------
int vtkOGSDepthProfile::RequestInformation( vtkInformation *vtkNotUsed(request),
	vtkInformationVector **inputVector, vtkInformationVector *outputVector) {

	// get the info objects
	vtkInformation *inInfo     = inputVector[0]->GetInformationObject(0);
	vtkInformation *sourceInfo = inputVector[1]->GetInformationObject(0);
	vtkInformation *outInfo    = outputVector->GetInformationObject(0);

	outInfo->CopyEntry(sourceInfo,vtkStreamingDemandDrivenPipeline::TIME_STEPS());
	outInfo->CopyEntry(sourceInfo,vtkStreamingDemandDrivenPipeline::TIME_RANGE());

	outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(),
				inInfo->Get(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT()),6);

	return 1;
}

//----------------------------------------------------------------------------
int vtkOGSDepthProfile::RequestUpdateExtent(vtkInformation *vtkNotUsed(request),
	vtkInformationVector **inputVector,vtkInformationVector *outputVector) {

	// get the info objects
	vtkInformation *inInfo     = inputVector[0]->GetInformationObject(0);
	vtkInformation *sourceInfo = inputVector[1]->GetInformationObject(0);
	vtkInformation *outInfo    = outputVector->GetInformationObject(0);

	int usePiece = 0;

	// What ever happened to CopyUpdateExtent in vtkDataObject?
	// Copying both piece and extent could be bad.  Setting the piece
	// of a structured data set will affect the extent.
	vtkDataObject* output = outInfo->Get(vtkDataObject::DATA_OBJECT());
	if (output &&
		(!strcmp(output->GetClassName(), "vtkUnstructuredGrid") ||
			!strcmp(output->GetClassName(), "vtkPolyData")))
		usePiece = 1;

	inInfo->Set(vtkStreamingDemandDrivenPipeline::EXACT_EXTENT(), 1);

	sourceInfo->Remove(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT());
	if (sourceInfo->Has(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT()))
		sourceInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(),
			sourceInfo->Get(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT()), 6);

	// SpatialMatch does not exist in this implementation
	sourceInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER(), 0);
	sourceInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES(), 1);
	sourceInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS(), 0);

	if (usePiece) {
		inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER(),
			outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER()));
		inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES(),
			outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES()));
		inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS(),
			outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS()));
	} else {
		inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(),
			outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT()), 6);
	}

	return 1;
}

//----------------------------------------------------------------------------
void vtkOGSDepthProfile::Initialize(vtkDataSet* input,vtkDataSet* source, vtkDataSet* output) {
	
	// Build Point List
	delete this->PointList;
	this->PointList = new vtkDataSetAttributes::FieldList(1);
	this->PointList->InitializeFieldList(source->GetPointData());

	// Build Cell List
	delete this->CellList;
	this->CellList = new vtkDataSetAttributes::FieldList(1);
	this->CellList->InitializeFieldList(source->GetCellData());

	vtkIdType numPts = input->GetNumberOfPoints();

	// Allocate storage for output PointData
	// All input PD is passed to output as PD. Those arrays in input CD that are
	// not present in output PD will be passed as output PD.
	output->GetPointData()->InterpolateAllocate((*this->PointList), numPts, numPts);
	// We assume we either have points or cells
	output->GetPointData()->InterpolateAllocate((*this->CellList), numPts, numPts);

	// Initialize Output Arrays
	int nArr = output->GetPointData()->GetNumberOfArrays();
	for (int arrId = 0; arrId < nArr; arrId++) {
		output->GetPointData()->GetArray(arrId)->SetNumberOfTuples(numPts);
		output->GetPointData()->GetArray(arrId)->Fill(0);
	}
}

void vtkOGSDepthProfile::Interpolate(vtkDataSet *input, vtkDataSet *source, vtkDataSet *output) {

	// Preallocate weights
	double *weights;
	weights = new double[source->GetMaxCellSize()];

	// Create the cell locator object
	vtkCellLocator *cellLocator = vtkCellLocator::New();
	cellLocator->SetDataSet(source); cellLocator->BuildLocator();

	// Recover the number of points in the input
	int npoints = input->GetNumberOfPoints();

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
			output->GetPointData()->InterpolatePoint(
					(*this->PointList), source->GetPointData(), 0, ii, cell->PointIds, weights
				);
			// Interpolate cell array data
			for (int arrId=0; arrId < source->GetCellData()->GetNumberOfArrays(); arrId++) {
				output->GetPointData()->CopyTuple(
					source->GetCellData()->GetArray(arrId),
					output->GetPointData()->GetArray(arrId),
					cellId,ii);
			}
		}
	}
	delete [] weights;
	cellLocator->Delete();
}

//----------------------------------------------------------------------------
void vtkOGSDepthProfile::SetSourceConnection(vtkAlgorithmOutput* algOutput) {
	this->SetInputConnection(1, algOutput);
}

void vtkOGSDepthProfile::SetSourceData(vtkDataObject *input) {
	this->SetInputData(1, input);
}

vtkDataObject *vtkOGSDepthProfile::GetSource() {
	if (this->GetNumberOfInputConnections(1) < 1)
		return nullptr;

	return this->GetExecutive()->GetInputData(1, 0);
}

//----------------------------------------------------------------------------
