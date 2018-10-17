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
#include "vtkInformation.h"
#include "vtkInformationVector.h"

#include "vtkObjectFactory.h"

vtkStandardNewMacro(vtkOGSTimeStatistics);


//----------------------------------------------------------------------------
vtkOGSTimeStatistics::vtkOGSTimeStatistics() {
}

//----------------------------------------------------------------------------
vtkOGSTimeStatistics::~vtkOGSTimeStatistics() {
}

//----------------------------------------------------------------------------
int vtkOGSTimeStatistics::RequestData(vtkInformation *vtkNotUsed(request), 
	vtkInformationVector **inputVector, vtkInformationVector *outputVector) {

	// get the info objects
	vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
	vtkInformation *outInfo = outputVector->GetInformationObject(0);

	// get the input and output
	vtkDataSet *input = vtkDataSet::SafeDownCast(
		inInfo->Get(vtkDataObject::DATA_OBJECT()));
	vtkDataSet *output = vtkDataSet::SafeDownCast(
		outInfo->Get(vtkDataObject::DATA_OBJECT()));

	output->CopyStructure(input);



	// Update progress and leave
	this->UpdateProgress(1.);
	return 1;
}

//----------------------------------------------------------------------------
int vtkOGSTimeStatistics::RequestInformation(vtkInformation* vtkNotUsed(request),
  vtkInformationVector** vtkNotUsed(inputVector), vtkInformationVector* outputVector) {


	return 1;
}
