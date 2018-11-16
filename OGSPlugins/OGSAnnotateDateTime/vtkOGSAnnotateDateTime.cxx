/*=========================================================================

  Program:   OGSAnnotateDateTime
  Module:    vtkOGSAnnotateDateTime.cxx

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkOGSAnnotateDateTime.h"

#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkDataSet.h"
#include "vtkFieldData.h"
#include "vtkStringArray.h"

#include "vtkObjectFactory.h"

#include <ctime>

vtkStandardNewMacro(vtkOGSAnnotateDateTime);

//----------------------------------------------------------------------------
vtkOGSAnnotateDateTime::vtkOGSAnnotateDateTime()
{
	this->TimeFormat = NULL;
}

//----------------------------------------------------------------------------
vtkOGSAnnotateDateTime::~vtkOGSAnnotateDateTime()
{
	this->SetTimeFormat(0);
}

//----------------------------------------------------------------------------
int vtkOGSAnnotateDateTime::RequestData(vtkInformation* request,
                                        vtkInformationVector** inputVector,
                                        vtkInformationVector* outputVector)
{

	// Recover the Date string
	vtkInformation *inputInfo = inputVector[0]->GetInformationObject(0);

	vtkDataSet *input = vtkDataSet::SafeDownCast(
		inputInfo->Get(vtkDataObject::DATA_OBJECT()));

	// Recover datevec
	vtkStringArray *vtkdate = vtkStringArray::SafeDownCast(
		input->GetFieldData()->GetAbstractArray("Date"));

	// If successful, modify the entry by parsing the string
	// using the time library
	if (vtkdate) {
		struct tm tm;
		char buff[256];
		strptime(vtkdate->GetValue(0).c_str(),"%Y%m%d-%H:%M:%S",&tm);
		strftime(buff,256,this->TimeFormat,&tm);
		this->Superclass::SetFormat(buff);
	}

	// Run RequestData from PythonAnnotationFilter
	return this->Superclass::RequestData(request,inputVector,outputVector);
}
