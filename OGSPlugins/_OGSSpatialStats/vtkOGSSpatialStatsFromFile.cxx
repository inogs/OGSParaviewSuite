/*=========================================================================

  Program:   OGSSpatialStats
  Module:    vtkOGSSpatialStatsFromFile.cxx

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkFloatArray.h"
#include "vtkCellData.h"
#include "vtkDataArraySelection.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkRectilinearGrid.h"

#include "vtkOGSSpatialStatsFromFile.h"

#include "vtkObjectFactory.h"

#define INDEX(ii,jj,kk,nx,ny) ( (nx-1)*(ny-1)*(kk) + (nx-1)*(jj) + (ii) )
#define GETVTKVAL1(vtkarray,ii,jj,kk,nx,ny)     ( (vtkarray)->GetTuple1((nx-1)*(ny-1)*(kk) + (nx-1)*(jj) + (ii)) )
#define GETVTKVAL3(vtkarray,ii,jj,kk,nx,ny)     ( (vtkarray)->GetTuple3((nx-1)*(ny-1)*(kk) + (nx-1)*(jj) + (ii)) )
#define SETVTKVAL1(vtkarray,val,ii,jj,kk,nx,ny) ( (vtkarray)->SetTuple1((nx-1)*(ny-1)*(kk) + (nx-1)*(jj) + (ii),(val)) )

vtkStandardNewMacro(vtkOGSSpatialStatsFromFile);


//----------------------------------------------------------------------------
vtkOGSSpatialStatsFromFile::vtkOGSSpatialStatsFromFile()
{
}

//----------------------------------------------------------------------------
vtkOGSSpatialStatsFromFile::~vtkOGSSpatialStatsFromFile()
{
}

//----------------------------------------------------------------------------
int vtkOGSSpatialStatsFromFile::RequestData(vtkInformation *vtkNotUsed(request),
	vtkInformationVector **inputVector,vtkInformationVector *outputVector) {
	// Get the info objects
	vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
	vtkInformation *outInfo = outputVector->GetInformationObject(0);

	// Get the input and output
	vtkRectilinearGrid *input = vtkRectilinearGrid::SafeDownCast(
		inInfo->Get(vtkDataObject::DATA_OBJECT()));
	vtkRectilinearGrid *output = vtkRectilinearGrid::SafeDownCast(
		outInfo->Get(vtkDataObject::DATA_OBJECT()));

	this->UpdateProgress(0.);

	
	
	// Copy the input grid
	this->UpdateProgress(1.);
	output->ShallowCopy(input);

	return 1;
}