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
#include "vtkStringArray.h"
#include "vtkCellData.h"
#include "vtkFieldData.h"
#include "vtkDataArraySelection.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkRectilinearGrid.h"

#include "vtkOGSSpatialStatsFromFile.h"

#include "vtkObjectFactory.h"

namespace NetCDF
{
// Include NetCDF functions
#include "vtknetcdf/include/netcdf.h"
}

#define INDEX(ii,jj,kk,nx,ny) ( (nx-1)*(ny-1)*(kk) + (nx-1)*(jj) + (ii) )
#define GETVTKVAL1(vtkarray,ii,jj,kk,nx,ny)     ( (vtkarray)->GetTuple1((nx-1)*(ny-1)*(kk) + (nx-1)*(jj) + (ii)) )
#define GETVTKVAL3(vtkarray,ii,jj,kk,nx,ny)     ( (vtkarray)->GetTuple3((nx-1)*(ny-1)*(kk) + (nx-1)*(jj) + (ii)) )
#define SETVTKVAL1(vtkarray,val,ii,jj,kk,nx,ny) ( (vtkarray)->SetTuple1((nx-1)*(ny-1)*(kk) + (nx-1)*(jj) + (ii),(val)) )

vtkStandardNewMacro(vtkOGSSpatialStatsFromFile);

//----------------------------------------------------------------------------
double *readNetCDF(const char *fname, const char *varname, int nb, int nc, int nz, int ns) {
	/*
		Read data from NetCDF files. Uses ParaView's built in NetCDF library.
	*/
	int fid, varid;
	double *out;

	// Open file for reading
	int retval = NetCDF::nc_open(fname,NC_NOWRITE,&fid);
	if ( retval != NC_NOERR )
		return NULL;
	// Get the variable id based on its name
	if ( NetCDF::nc_inq_varid(fid,varname,&varid) != NC_NOERR )
		return NULL;
	// Read the data
	out = (double*)malloc(nb*nc*nz*ns*sizeof(double));
	NetCDF::nc_get_var_double(fid,varid,out);
	// Close the file
	NetCDF::nc_close(fid);
	// Eliminate the missing variables
//	for (int ii=0;ii<nx*ny*nz;ii++)
//		if(out[ii] > MAXVAL) out[ii] = 0.;

	// Return
	return out;
}

//----------------------------------------------------------------------------
void addStat(vtkDataArraySelection *Array) {
	Array->AddArray("Mean");
	Array->AddArray("std");
	Array->AddArray("min");
	Array->AddArray("p05");
	Array->AddArray("p25");
	Array->AddArray("p50");
	Array->AddArray("p75");
	Array->AddArray("p95");
	Array->AddArray("max");
}

//----------------------------------------------------------------------------
vtkOGSSpatialStatsFromFile::vtkOGSSpatialStatsFromFile(){
	this->StatDataArraySelection = vtkDataArraySelection::New();
	addStat(this->StatDataArraySelection);

	this->FolderName  = NULL;
	this->bmask_field = NULL;
	this->cmask_field = NULL;
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

	// Loop the number of variables
	// For each variable, we will see if a prost processing exists and
	// then we will loop the mesh and create the statistics.
	int narrays = input->GetCellData()->GetNumberOfArrays();
	int nstat   = this->GetNumberOfStatArrays();
	for (int varId = 0; varId < narrays; varId++) {
		// Recover the array and the array name
		vtkFloatArray *vtkVarArray = vtkFloatArray::SafeDownCast(
			input->GetCellData()->GetArray(varId));
		char *array_name = vtkVarArray->GetName();

		// Do not work with the basins and coasts mask
		if (std::string(basins_mask->GetName()) == std::string(array_name)) continue;
		if (std::string(coasts_mask->GetName()) == std::string(array_name)) continue;

		// At this point, we can try to load the stat profile
		double *stat_profile = readNetCDF(filename,    // Filename to read
										  array_name,  // Array to obtain
										  nbasins,     // Number of sub basins
										  ncoasts,     // Number of coasts
										  nz,		   // Number of z coords
										  nstat        // Number of statistics
										 );
		// If file cannot be read or variable does not exist
		if (stat_profile == NULL) {
			vtkWarningMacro("File <"<<filename<<"> or variable <"
				<<array_name<<"> cannot be read!");
			continue;
		}

		// At this point we do have the statistics per basin, coast and depth level
		// of a single variable. We must initialize a vector of vtkFloatArrays where
		// the statistics will be stored during the loop.

printf("%s\n",array_name);
	}
	
	
	// Copy the input grid
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
