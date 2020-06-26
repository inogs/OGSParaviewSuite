/*=========================================================================

  Program:   OGSTimeCommons
  Module:    vtkOGSTimeCommons.cxx

  Copyright (c) 2020 Arnau Miro, OGS
  All rights reserved.

	 This software is distributed WITHOUT ANY WARRANTY; without even
	 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
	 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

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

#include <vector>
#include <string>
#include <ctime>
#include <chrono>

#include "TimeList.h"

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

//----------------------------------------------------------------------------
void RecoverMasterFileName(std::string &fname, vtkDataSet *input) {

	// Recover the master file name from the metadata array
	// Return whether we need to stop executing or not

	vtkStringArray *vtkmetadata = vtkStringArray::SafeDownCast(
		input->GetFieldData()->GetAbstractArray("Metadata"));

	// If successful, recover the file name
	if (vtkmetadata)
		fname = vtkmetadata->GetValue(4);
	else
		fname = std::string("NotFound");
}