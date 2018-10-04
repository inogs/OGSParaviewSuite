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
#include "vtkCharArray.h"
#include "vtkGenericCell.h"
#include "vtkPointData.h"
#include "vtkPointSet.h"
#include "vtkDataSet.h"
#include "vtkImageData.h"
#include "vtkExecutive.h"
#include "vtkSmartPointer.h"
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

	this->MaskPoints = nullptr;
	this->PointList  = nullptr;
	this->CellList   = nullptr;

	this->CellLocatorPrototype = nullptr;

	this->ValidPointMaskArrayName = nullptr;
	this->SetValidPointMaskArrayName("vtkValidPointMask");
}

//----------------------------------------------------------------------------
vtkOGSDepthProfile::~vtkOGSDepthProfile() {
	delete this->PointList;
	delete this->CellList;

	this->vtkOGSDepthProfile::SetValidPointMaskArrayName(nullptr);
	this->vtkOGSDepthProfile::SetCellLocatorPrototype(nullptr);

	if (this->MaskPoints)
		this->MaskPoints->Delete();
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

	// If there is data to interpolate, begin the interpolation
	if (source) {
		this->BuildFieldList(source);
		this->InitializeForProbing(input, output);
		this->Interpolate(input, source, output);
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

	// A variation of the bug fix from John Biddiscombe.
	// Make sure that the scalar type and number of components
	// are propagated from the source not the input.
	if (vtkImageData::HasScalarType(sourceInfo))
		vtkImageData::SetScalarType(vtkImageData::GetScalarType(sourceInfo),outInfo);
	
	if (vtkImageData::HasNumberOfScalarComponents(sourceInfo))
		vtkImageData::SetNumberOfScalarComponents(vtkImageData::GetNumberOfScalarComponents(sourceInfo),outInfo);

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
void vtkOGSDepthProfile::BuildFieldList(vtkDataSet* source) {
  delete this->PointList;
  delete this->CellList;

  this->PointList = new vtkDataSetAttributes::FieldList(1);
  this->PointList->InitializeFieldList(source->GetPointData());

  this->CellList = new vtkDataSetAttributes::FieldList(1);
  this->CellList->InitializeFieldList(source->GetCellData());
}

void vtkOGSDepthProfile::InitializeForProbing(vtkDataSet* input,vtkDataSet* output) {

	vtkIdType numPts = input->GetNumberOfPoints();

	// if this is repeatedly called by the pipeline for a composite mesh,
	// you need a new array for each block
	// (that is you need to reinitialize the object)
	if (this->MaskPoints)
		this->MaskPoints->Delete();

	this->MaskPoints = vtkCharArray::New();
	this->MaskPoints->SetNumberOfComponents(1);
	this->MaskPoints->SetNumberOfTuples(numPts);
	this->MaskPoints->FillValue(0);
	this->MaskPoints->SetName(this->ValidPointMaskArrayName ? 
		this->ValidPointMaskArrayName : "vtkValidPointMask");

	// Allocate storage for output PointData
	// All input PD is passed to output as PD. Those arrays in input CD that are
	// not present in output PD will be passed as output PD.
	vtkPointData* outPD = output->GetPointData();
	outPD->InterpolateAllocate((*this->PointList), numPts, numPts);
	// We assume we either have points or cells
	outPD->InterpolateAllocate((*this->CellList), numPts, numPts);

	this->InitializeOutputArrays(outPD, numPts);
	outPD->AddArray(this->MaskPoints);
}

void vtkOGSDepthProfile::InitializeOutputArrays(vtkPointData *outPD, vtkIdType numPts) {
	for (int i = 0; i < outPD->GetNumberOfArrays(); ++i) {
		vtkDataArray* da = outPD->GetArray(i);
		if (da) {
			da->SetNumberOfTuples(numPts);
			da->Fill(0);
		}
	}
}

void vtkOGSDepthProfile::Interpolate(vtkDataSet *input, vtkDataSet *source, vtkDataSet *output) {

	int subId;
	double xyz[3], pcoords[3], closestPoint[3], *weights;

	double tol2 = 400;

	// Preallocate weights
	weights = new double[source->GetMaxCellSize()];

	// Create the cell locator object
	vtkSmartPointer<vtkCellLocator> cellLocator = vtkSmartPointer<vtkCellLocator>::New();
	cellLocator->SetDataSet(source);
	cellLocator->BuildLocator();

	char* maskArray = this->MaskPoints->GetPointer(0);

	// Recover the number of points in the input
	int npoints = input->GetNumberOfPoints();

	// Loop the number of points
	vtkNew<vtkGenericCell> gcell;
	int abort = 0;
	for (int ii = 0; ii < npoints && !abort; ii++) {
		// Update progress
		this->UpdateProgress((double)(ii)/(double)(npoints));
		abort = GetAbortExecute();

		// skip points which have already been probed with success.
		// This is helpful for multiblock dataset probing.
		if (maskArray[ii] == static_cast<char>(1))
			continue;

		// Get the xyz coordinate of the point in the input dataset
		input->GetPoint(ii, xyz);

		// Find the cell that contains xyz and get it
		vtkIdType cellId;
		double dist2;
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
			output->GetPointData()->InterpolatePoint(
					(*this->CellList), source->GetCellData(), 0, ii, cell->PointIds, weights
				);			
			// Update mask array
			maskArray[ii] = static_cast<char>(1);
		}
	}

	this->MaskPoints->Modified();
	delete [] weights;
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
