/*=========================================================================

  Program:   OGSSpatialStats
  Module:    vtkOGSSpatialStatsFromFile.cxx

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkCell.h"
#include "vtkPoints.h"
#include "vtkFloatArray.h"
#include "vtkStringArray.h"
#include "vtkCellData.h"
#include "vtkFieldData.h"
#include "vtkDataArraySelection.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#include "vtkOGSSpatialStatsFromFile.h"

#include "vtkObjectFactory.h"

#include <string>
#include <algorithm>
#include <map>


namespace NetCDF
{
// Include NetCDF IO
#include "../_utils/netcdfio.cpp"
}

namespace VTK
{
// Include the VTK Operations
#include "../_utils/vtkOperations.cpp"
}


vtkStandardNewMacro(vtkOGSSpatialStatsFromFile);

//----------------------------------------------------------------------------
vtkOGSSpatialStatsFromFile::vtkOGSSpatialStatsFromFile(){
	this->StatDataArraySelection = vtkDataArraySelection::New();
	this->StatDataArraySelection->AddArray("Mean");
	this->StatDataArraySelection->AddArray("std");
	this->StatDataArraySelection->AddArray("min");
	this->StatDataArraySelection->AddArray("p05");
	this->StatDataArraySelection->AddArray("p25");
	this->StatDataArraySelection->AddArray("p50");
	this->StatDataArraySelection->AddArray("p75");
	this->StatDataArraySelection->AddArray("p95");
	this->StatDataArraySelection->AddArray("max");

	this->FolderName  = NULL;
	this->bmask_field = NULL;
	this->cmask_field = NULL;

	this->per_coast = 0;
}

//----------------------------------------------------------------------------
vtkOGSSpatialStatsFromFile::~vtkOGSSpatialStatsFromFile() {
	this->StatDataArraySelection->Delete();

	this->SetFolderName(0);
	this->Setbmask_field(0);
	this->Setcmask_field(0);
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

	// We just want to copy the mesh, not the variables
	output->CopyStructure(input); 

	this->UpdateProgress(0.);

	// First we need to recover the datetime inside the input
	vtkStringArray *strf = vtkStringArray::SafeDownCast( 
		input->GetFieldData()->GetAbstractArray("Date") );
	const char *datetime = strf->GetValue(0);
	// Construct the file name to open along with the path
	char filename[256];
	sprintf(filename,"%s/ave.%s.stat_profiles.nc",this->FolderName,datetime);

	// Next, recover the dimensions of the rectilinear grid
	vtkFloatArray *vtkxcoord = vtkFloatArray::SafeDownCast(
		input->GetXCoordinates());
	int nx = vtkxcoord->GetNumberOfTuples();
	vtkFloatArray *vtkycoord = vtkFloatArray::SafeDownCast(
		input->GetYCoordinates());
	int ny = vtkycoord->GetNumberOfTuples();
	vtkFloatArray *vtkzcoord = vtkFloatArray::SafeDownCast(
		input->GetZCoordinates());
	int nz = vtkzcoord->GetNumberOfTuples();

	// Also recover the basins and coasts mask
	// and add them to the output
	vtkFloatArray *basins_mask = vtkFloatArray::SafeDownCast(
		input->GetCellData()->GetArray(this->bmask_field));
	output->GetCellData()->AddArray(basins_mask);
	int nbasins = 21;
	vtkFloatArray *coasts_mask = vtkFloatArray::SafeDownCast(
		input->GetCellData()->GetArray(this->cmask_field));
	output->GetCellData()->AddArray(coasts_mask);
	int ncoasts = 3;

	// Also copy e1t and e2t
	vtkFloatArray *vtke1t = vtkFloatArray::SafeDownCast(
		input->GetCellData()->GetArray("e1t"));
	output->GetCellData()->AddArray(vtke1t);
	vtkFloatArray *vtke2t = vtkFloatArray::SafeDownCast(
		input->GetCellData()->GetArray("e2t"));
	output->GetCellData()->AddArray(vtke2t);

	// Loop the number of variables
	// For each variable, we will see if a post processing exists and
	// then we will loop the mesh and create the statistics.
	int narrays = input->GetCellData()->GetNumberOfArrays();
	int nstat   = this->GetNumberOfStatArrays();
	for (int varId = 0; varId < narrays; varId++) {
		// Recover the array and the array name
		vtkFloatArray *vtkVarArray = vtkFloatArray::SafeDownCast(
			input->GetCellData()->GetArray(varId));
		char *array_name = vtkVarArray->GetName();

		// Do not work with the basins, coasts mask, e1t or e2t
		if (std::string(basins_mask->GetName()) == std::string(array_name)) continue;
		if (std::string(coasts_mask->GetName()) == std::string(array_name)) continue;
		if (std::string(vtke1t->GetName())      == std::string(array_name)) continue;
		if (std::string(vtke2t->GetName())      == std::string(array_name)) continue;

		// At this point, we can try to load the stat profile
		double *stat_profile = NetCDF::readNetCDF(filename, array_name, nbasins*ncoasts*nz*nstat);
		// If file cannot be read or variable does not exist
		if (stat_profile == NULL) {
			vtkWarningMacro("File <"<<filename<<"> or variable <"
				<<array_name<<"> cannot be read!");
			continue;
		}

		// At this point we do have the statistics per basin, coast and depth level
		// of a single variable. We must initialize a vector of vtkFloatArrays where
		// the statistics will be stored during the loop.
		std::map<std::string, vtkFloatArray*> mapStatArray;

		for (int statId = 0; statId < nstat; statId++) {
			const char *statName = this->GetStatArrayName(statId);
			// Skip those arrays who have not been enabled
			if (!this->GetStatArrayStatus(statName))
				continue;
			// Statistical variable name
			char statVarName[256];
			sprintf(statVarName,"%s, %s",array_name,statName);
			// Define a new vtkFloatArray
			vtkFloatArray *vtkStatVar = VTK::createVTKscaf(statVarName, nx, ny, nz, NULL);
			// Store the array in the map
			mapStatArray.insert(std::make_pair(std::string(statName),vtkStatVar));
		}

		// We are now ready to loop the mesh and set the variables accordingly
		for (int kk = 0; kk < nz; kk++) {
			for (int jj = 0; jj < ny; jj++) {
				for (int ii = 0; ii < nx; ii++) {
					// Which basin and which coast are we?
					int basinId = (int)( basins_mask->GetTuple1(VTKIND(ii,jj,kk,nx,ny)) ) - 1;
					int coastId = this->per_coast ? 
						(int)( coasts_mask->GetTuple1(VTKIND(ii,jj,kk,nx,ny)) ) - 1 : 2;
					// Loop all the requested statistics using an iterator
					std::map<std::string,vtkFloatArray*>::iterator iter;
					for (iter = mapStatArray.begin(); iter != mapStatArray.end(); iter++) {
						// Which statistic are we computing?
						int statId = this->GetStatArrayIndex(iter->first.c_str());
						// Recover value from the stat profile
						double value = stat_profile[STTIND(basinId,coastId,kk,statId,nstat,nz,ncoasts)];
						// Set the value
						iter->second->SetTuple1(VTKIND(ii,jj,kk,nx,ny),value);
					}
				}
			}
		}

		// Now that we computed the arrays, we can set them in the output
		// and deallocate memory
		std::map<std::string,vtkFloatArray*>::iterator iter;
		for (iter = mapStatArray.begin(); iter != mapStatArray.end(); iter++) {
			output->GetCellData()->AddArray(iter->second);
			iter->second->Delete();
		}

		this->UpdateProgress(0.+1./(double)(narrays)*(double)(varId));
	}
	
	// Update progress and exit successfully
	this->UpdateProgress(1.);
	return 1;
}

//----------------------------------------------------------------------------
void vtkOGSSpatialStatsFromFile::DisableAllStatArrays()
{
	this->StatDataArraySelection->DisableAllArrays();
}

void vtkOGSSpatialStatsFromFile::EnableAllStatArrays()
{
	this->StatDataArraySelection->EnableAllArrays();
}

int vtkOGSSpatialStatsFromFile::GetNumberOfStatArrays()
{
	return this->StatDataArraySelection->GetNumberOfArrays();
}

const char* vtkOGSSpatialStatsFromFile::GetStatArrayName(int index)
{
	if (index >= (int)this->GetNumberOfStatArrays() || index < 0)
		return NULL;
	else
		return this->StatDataArraySelection->GetArrayName(index);
}

int vtkOGSSpatialStatsFromFile::GetStatArrayIndex(const char* name)
{
	return this->StatDataArraySelection->GetArrayIndex(name);
}

int vtkOGSSpatialStatsFromFile::GetStatArrayStatus(const char* name)
{
	return this->StatDataArraySelection->ArrayIsEnabled(name);
}

void vtkOGSSpatialStatsFromFile::SetStatArrayStatus(const char* name, int status)
{
	if (status)
		this->StatDataArraySelection->EnableArray(name);
	else
		this->StatDataArraySelection->DisableArray(name);

	this->Modified();
}
