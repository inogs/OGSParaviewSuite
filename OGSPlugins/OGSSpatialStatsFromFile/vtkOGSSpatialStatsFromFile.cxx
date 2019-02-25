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
#include "vtkTypeUInt8Array.h"
#include "vtkFloatArray.h"
#include "vtkDoubleArray.h"
#include "vtkStringArray.h"
#include "vtkCellData.h"
#include "vtkFieldData.h"
#include "vtkDataArray.h"
#include "vtkDataArraySelection.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#include "vtkOGSSpatialStatsFromFile.h"

#include "vtkObjectFactory.h"

#include <cstdint>
#include <string>

vtkStandardNewMacro(vtkOGSSpatialStatsFromFile);

//----------------------------------------------------------------------------
/*
	Macro to set the array precision 
*/
#define FLDARRAY double
#define VTKARRAY vtkDoubleArray

#define STTIND(bId,cId,kk,sId,ns,nz,nc) ( (ns)*(nz)*(nc)*(bId) + (ns)*(nz)*(cId) + (ns)*(kk) + (sId) )

#include "../_utils/field.h"
#include "../_utils/vtkFields.hpp"
#include "../_utils/netcdfio.hpp"

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

	this->SetFolderName(NULL);
	this->Setbmask_field(NULL);
	this->Setcmask_field(NULL);
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

	// Copy Metadata array
	vtkStringArray *vtkmetadata = vtkStringArray::SafeDownCast(
		input->GetFieldData()->GetAbstractArray("Metadata"));
	output->GetFieldData()->AddArray(vtkmetadata);

	this->UpdateProgress(0.);

	// First we need to recover the date from the metadata
	std::string filename = std::string(this->FolderName) + 
		                   std::string("/ave.") + 
		                   vtkmetadata->GetValue(0) + 
		                   std::string(".stat_profiles.nc");

	// Next, recover the dimensions of the rectilinear grid
	int nx = input->GetXCoordinates()->GetNumberOfTuples();
	int ny = input->GetYCoordinates()->GetNumberOfTuples();
	int nz = input->GetZCoordinates()->GetNumberOfTuples();

	// Recover the basins and coasts mask
	// and add them to the output
	vtkTypeUInt8Array *vtkMask;
	int nbasins = 21, ncoasts = 3;
	field::Field<uint8_t> bmask, cmask;
	// Basins mask
	vtkMask = vtkTypeUInt8Array::SafeDownCast( 
		input->GetCellData()->GetArray(this->bmask_field) );
	output->GetCellData()->AddArray(vtkMask);
	bmask = VTK::createFieldfromVTK<vtkTypeUInt8Array,uint8_t>(vtkMask);
	// Coast mask
	vtkMask = vtkTypeUInt8Array::SafeDownCast( 
		input->GetCellData()->GetArray(this->cmask_field) );
	output->GetCellData()->AddArray(vtkMask);
	cmask = VTK::createFieldfromVTK<vtkTypeUInt8Array,uint8_t>(vtkMask);

	// Copy e1, e2 and e3
	VTKARRAY *vtkArray;
	vtkArray = VTKARRAY::SafeDownCast( 
		input->GetCellData()->GetArray("e1") );
	output->GetCellData()->AddArray(vtkArray);
	vtkArray = VTKARRAY::SafeDownCast( 
		input->GetCellData()->GetArray("e2") );
	output->GetCellData()->AddArray(vtkArray);
		vtkArray = VTKARRAY::SafeDownCast( 
		input->GetCellData()->GetArray("e3") );
	output->GetCellData()->AddArray(vtkArray);

	// Loop the number of variables
	// For each variable, we will see if a post processing exists and
	// then we will loop the mesh and create the statistics.
	int narrays = input->GetCellData()->GetNumberOfArrays();
	int nstat   = this->GetNumberOfStatArrays();

	for (int varId = 0; varId < narrays; varId++) {
		// Recover the array and the array name
		vtkDataArray *vtkDArray;
		vtkDArray = input->GetCellData()->GetArray(varId);
		std::string arrName = vtkDArray->GetName();

		// Do not work with the basins, coasts mask, e1, e2 or e3
		if (std::string(this->bmask_field) == arrName) continue;
		if (std::string(this->cmask_field) == arrName) continue;
		if (std::string("e1")              == arrName) continue;
		if (std::string("e2")              == arrName) continue;
		if (std::string("e3")              == arrName) continue;

		// Recover Array values
		field::Field<FLDARRAY> array;
		vtkArray = VTKARRAY::SafeDownCast( vtkDArray );
		array = VTK::createFieldfromVTK<VTKARRAY,FLDARRAY>(vtkArray);

		// At this point, we can try to load the stat profile
		field::Field<FLDARRAY> statProfile(nbasins*ncoasts*(nz-1)*nstat,1);
		
		if ( NetCDF::readNetCDF2F(filename.c_str(),arrName.c_str(),statProfile) != NETCDF_OK ) {
			// If file cannot be read or variable does not exist
			vtkWarningMacro("File <"<<filename.c_str()<<"> or variable <"
				<<arrName.c_str()<<"> cannot be read!");
			continue;
		}

		// At this point we do have the statistics per basin, coast and depth level
		// of a single variable. Loop the number of statistics.
		// Parallelization MPI
		for (int statId = 0; statId < nstat; ++statId) {
			// Recover stat name
			std::string statName = this->GetStatArrayName(statId);
			
			// Skip those arrays who have not been enabled
			if (!this->GetStatArrayStatus(statName.c_str()))
				continue;
			
			// Generate variable name
			statName = arrName + std::string(", ") + statName;
			
			// Loop the mesh and generate the array
			field::Field<FLDARRAY> statArray(array.get_n(),array.get_m());
			// Parallelization OPENMP
			for (int kk = 0; kk < nz-1; kk++) {
				for (int jj = 0; jj < ny-1; jj++) {
					for (int ii = 0; ii < nx-1; ii++) {
						// Position acording x,y,z
						int pos = CLLIND(ii,jj,kk,nx,ny);
						// In which basin are we? (we need to loop the basins and find which is true)
						int  bId = 0; 
						bool isbasin = false;
						for (bId = 0; bId < bmask.get_m(); ++bId) {
							if (bmask[pos][bId]) { isbasin = true; break; }
						}
//						int bId = (int)( bmask[pos][0] ) - 1;
						// In which coast are we?
						int cId = this->per_coast ? cmask[pos][0] - 1 : 2;
						// Set the value (generally cId < 0 when bId < 0)
						statArray[pos][0] = (isbasin) ? statProfile[STTIND(bId,cId,kk,statId,nstat,nz-1,ncoasts)][0] : 0.;
					}
				}
			}

			// Set variable in the mesh and deallocate if needed
			VTKARRAY *vtkstatArray;
			vtkstatArray = VTK::createVTKfromField<VTKARRAY,FLDARRAY>(statName.c_str(),statArray);
			output->GetCellData()->AddArray(vtkstatArray);
			vtkstatArray->Delete();
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
