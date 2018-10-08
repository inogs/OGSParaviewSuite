/*=========================================================================

  Program:   OGSReader
  Module:    vtkOGSReader.cxx

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include <cstdlib>
#include <cstdio>
#include <string>
#include <cstring>

#include <vtksys/SystemTools.hxx>

#include "vtkCallbackCommand.h"
#include "vtkCommand.h"
#include "vtkFloatArray.h"
#include "vtkStringArray.h"
#include "vtkCellData.h"
#include "vtkDataArraySelection.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkRectilinearGrid.h"
#include "vtkObjectFactory.h"

#include "vtkOGSReader.h"

#ifdef PARAVIEW_USE_MPI
#include "vtkMultiProcessController.h"
vtkCxxSetObjectMacro(vtkOGSReader, Controller, vtkMultiProcessController);
#endif

#define INDEX(array,ii,jj,kk,nx,ny) ( (array)[(nx)*(ny)*(kk) + (nx)*(jj) + (ii)] )

vtkStandardNewMacro(vtkOGSReader);

//----------------------------------------------------------------------------

namespace OGS
{
// Include OGS specific functions
#include "../_utils/OGSmesh.cpp"
#include "../_utils/OGSfile.cpp"
}

namespace NetCDF
{
// Include NetCDF IO
#include "../_utils/netcdfio.cpp"
}


//----------------------------------------------------------------------------
void createRectilinearGrid(int nLon, int nLat, int nLev,
	double *Lon2Meters, double* Lat2Meters, double *nav_lev, double dps,
	vtkRectilinearGrid *rgrid) {
	/*
		Create a vtkRectilinearGrid given the mesh dimensions.
	*/
	// Set vtkFloatArrays
	vtkFloatArray *daLon = vtkFloatArray::New();
	for (int ii = 0; ii < nLon; ii++) daLon->InsertNextValue(1.*Lon2Meters[ii]);

	vtkFloatArray *daLat = vtkFloatArray::New();
	for (int ii = 0; ii < nLat; ii++) daLat->InsertNextValue(1.*Lat2Meters[ii]);

	vtkFloatArray *daLev = vtkFloatArray::New();
	for (int ii=0; ii<nLev; ii++) daLev->InsertNextValue(-dps*nav_lev[ii]);

	// Set rectilinear grid
	rgrid->SetDimensions(nLon,nLat,nLev);
	rgrid->SetXCoordinates(daLon);
	rgrid->SetYCoordinates(daLat);
	rgrid->SetZCoordinates(daLev);
}

vtkFloatArray *createVTKscaf(const char *name,int nx, int ny, int nz, double *array) {
	/*
		Create a vtk scalar field given the name, dimensions and
		the array to be converted.
	*/
	// Create float array
	vtkFloatArray *vtkArray = vtkFloatArray::New();
	vtkArray->SetName(name);
	vtkArray->SetNumberOfComponents(1); // Scalar field
	vtkArray->SetNumberOfTuples(nx*ny*nz);

	// Preallocate vtkArray to zero
	vtkArray->Fill(0.);
//	for (int ii = 0; ii < nx*ny*nz; ii++) 
//		vtkArray->InsertNextValue(0.);

	// Fill the vtkArray with the values of the array
	for (int kk = 0; kk < nz-1; kk++) {
		for (int jj = 0; jj < ny-1; jj++) {
			for (int ii = 0; ii < nx-1; ii++) {
				vtkArray->SetTuple1(ii+jj*(nx-1)+kk*(nx-1)*(ny-1),
					INDEX(array,ii,jj,kk,nx,ny));
			}
		}
	}

	return vtkArray;
}

vtkFloatArray *createVTKvecf3(const char *name,int nx, int ny, int nz, double *a1, double *a2, double *a3) {
	/*
		Create a vtk vectorial field given the name, dimensions and
		the array to be converted.
	*/
	// Create float array
	vtkFloatArray *vtkArray = vtkFloatArray::New();
	vtkArray->SetName(name);
	vtkArray->SetNumberOfComponents(3); // Scalar field
	vtkArray->SetNumberOfTuples(nx*ny*nz);

	// Preallocate vtkArray to zero
	vtkArray->Fill(0.);
//	for (int ii = 0; ii < 3*nx*ny*nz; ii++) 
//		vtkArray->InsertNextValue(0.);

	// Fill the vtkArray with the values of the array
	for (int kk = 0; kk < nz-1; kk++) {
		for (int jj = 0; jj < ny-1; jj++) {
			for (int ii = 0; ii < nx-1; ii++) {
				vtkArray->SetTuple3(ii+jj*(nx-1)+kk*(nx-1)*(ny-1),
					INDEX(a1,ii,jj,kk,nx,ny),
					INDEX(a2,ii,jj,kk,nx,ny),
					INDEX(a3,ii,jj,kk,nx,ny));
			}
		}
	}

	return vtkArray;
}

vtkStringArray *createVTKstrf(const char *name,const char *data) {
	/*
		Create a vtk string array given the name and data to be
		converted.
	*/
	// Create string array
	vtkStringArray *vtkArray = vtkStringArray::New();
	vtkArray->SetName(name);
	vtkArray->SetNumberOfTuples(1);

	// Set the value
	vtkArray->SetValue(0,data);

	return vtkArray;
}

//----------------------------------------------------------------------------
vtkOGSReader::vtkOGSReader() {

	this->SetNumberOfInputPorts(0);
	this->SetNumberOfOutputPorts(1);

	this->SubBasinsMask = 1;
	this->CoastsMask    = 1;
	this->DepthScale    = 1000.;

	this->Mesh = (vtkRectilinearGrid*)vtkRectilinearGrid::New();

	this->NumberOfAvePhysFields     = 0;
	this->NumberOfAvePhysComponents = 0;
	this->NumberOfAveFreqFields     = 0;
	this->NumberOfAveFreqComponents = 0;

	this->AvePhysDataArraySelection = vtkDataArraySelection::New();
	this->AveFreqDataArraySelection = vtkDataArraySelection::New();

	#ifdef PARAVIEW_USE_MPI
		this->Controller = NULL;
		this->SetController(vtkMultiProcessController::GetGlobalController());
	#endif
}

//----------------------------------------------------------------------------
vtkOGSReader::~vtkOGSReader() {

	this->FileName = NULL;
	this->SetFileName(0);

	if (this->ave_phys.nvars > 0) free(this->ave_phys.vars);
	if (this->ave_freq.nvars > 0) free(this->ave_freq.vars);

	if (this->timeStepInfo.ntsteps > 0) {
		for (int ii = 0; ii < this->timeStepInfo.ntsteps; ii++)
			free(this->timeStepInfo.datetime[ii]);
		free(this->timeStepInfo.datetime);
	}

	this->AvePhysDataArraySelection->Delete();
	this->AveFreqDataArraySelection->Delete();

	if(this->Mesh) this->DeleteMesh();

	#ifdef PARAVIEW_USE_MPI
  		this->SetController(NULL);	
	#endif
}

//----------------------------------------------------------------------------
int vtkOGSReader::RequestData(vtkInformation* vtkNotUsed(request),
  vtkInformationVector** vtkNotUsed(inputVector), vtkInformationVector* outputVector) {


	#ifdef PARAVIEW_USE_MPI
	  	if (this->Controller) {
	  		if (this->Controller->GetNumberOfProcesses() > 1)
	  			vtkWarningMacro("OGSReader is not parallel-aware: all pvserver processes will read the entire file!");
	  	}
	#endif

	vtkDebugMacro("Opening file: " << this->FileName);

  	// get the info object
	vtkInformation* outInfo = outputVector->GetInformationObject(0);

	// get the output
	vtkRectilinearGrid *output = 
	vtkRectilinearGrid::SafeDownCast(outInfo->Get(vtkRectilinearGrid::DATA_OBJECT()));

	// Clear the mesh
	if (this->Mesh) this->DeleteMesh();
	this->Mesh = (vtkRectilinearGrid*)vtkRectilinearGrid::New();

	this->UpdateProgress(0.0);

	/* SET THE TIME STEPPING

		The number of time steps are read in the master file. If that number is
		greater than one, then use the ParaView built in functions to move
		through the time steps.

	*/
	double requestedTimeValue = 0.;
	if (outInfo->Has(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP())) {
	    // Get the requested time step. We only support requests of a single time
	    // step in this reader right now
	    requestedTimeValue = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());

	    output->GetInformation()->Set(vtkDataObject::DATA_TIME_STEP(), requestedTimeValue);
	}
	int ii_tstep = 0;
	if (this->timeStepInfo.ntsteps > 1)
		ii_tstep = (int)requestedTimeValue;

	/* READING THE MESH FILE 
		
		The mesh data (dims, Lon2Meters, Lat2Meters and nav_lev) are inside
		the .ogsmsh file containing the correct resolution of the meshmask for
		the current simulation.

		See OGSmesh for further details.

	*/
	int nLon, nLat, nLev;
	double *Lon2Meters, *Lat2Meters, *nav_lev, *basins_mask, *coast_mask;

	OGS::readOGSMesh(this->meshfile,
					 &nLon, &nLat, &nLev,
					 &Lon2Meters, &Lat2Meters, &nav_lev,
					 &basins_mask, &coast_mask
					 );

	this->UpdateProgress(0.05);

	/* CREATE RECTILINEAR GRID

		The rectilinear grid is created here from the mesh data that has
		been previously read. Basins mask and coast mask are added according
		to user selection.

	*/
	int ncells = nLon*nLat*nLev;
	createRectilinearGrid(nLon,nLat,nLev,Lon2Meters,Lat2Meters,nav_lev,this->DepthScale,this->Mesh);
	free(Lon2Meters); free(Lat2Meters); free(nav_lev);

	// Sub-basins mask
	if (this->SubBasinsMask) {
		vtkFloatArray *vtkarray = createVTKscaf("basins mask",nLon,nLat,nLev,basins_mask);
		this->Mesh->GetCellData()->AddArray(vtkarray);
		vtkarray->Delete();
	}
	free(basins_mask);
	
	// Coasts mask
	if (this->CoastsMask) {
		vtkFloatArray *vtkarray = createVTKscaf("coast mask" ,nLon,nLat,nLev,coast_mask);
		this->Mesh->GetCellData()->AddArray(vtkarray);
		vtkarray->Delete();
	}
	free(coast_mask);

	this->UpdateProgress(0.25);

	/* READING THE PHYSICAL VARIABLES

		Variables inside AVE_PHYS are read here. Velocity is outputed as
		a vector while the others are scalar arrays. User can select which
		variables to load by using a panel.

	*/
	for (int ii = 0; ii < this->ave_phys.nvars; ii++) {
		char *varname = this->ave_phys.vars[ii].name;
		char *cdfname = this->ave_phys.vars[ii].vname;
		char *varpath = OGS::writeOGSPath(this->ave_phys.vars[ii].path,
			this->timeStepInfo.datetime[ii_tstep],"*");
		// Test if the variable has been activated
		if (this->GetAvePhysArrayStatus(varname)) {
			// Velocity is treated separately
			if (std::string(varname) == "Velocity") {
				double *u = NetCDF::readNetCDF(varpath,"vozocrtx",nLon*nLat*nLev),
					   *v = NetCDF::readNetCDF(varpath,"vomecrty",nLon*nLat*nLev),
					   *w = NetCDF::readNetCDF(varpath,"vovecrtz",nLon*nLat*nLev);
				vtkFloatArray *vtkarray = createVTKvecf3(varname,nLon,nLat,nLev,u,v,w);
				this->Mesh->GetCellData()->AddArray(vtkarray);
				vtkarray->Delete(); free(u); free(v); free(w);
			}
			else {
				double *array = NetCDF::readNetCDF(varpath,cdfname,nLon*nLat*nLev);
				vtkFloatArray *vtkarray = createVTKscaf(varname,nLon,nLat,nLev,array);
				this->Mesh->GetCellData()->AddArray(vtkarray);
				vtkarray->Delete(); free(array);
			}
		}
		free(varpath);
		this->UpdateProgress(0.25 + ii*(0.375/ave_phys.nvars));		
	}

	/* READING THE BIOGEOCHEMICAL VARIABLES

		Variables inside AVE_FREQ are read here. User can select which
		variables to load by using a panel.

	*/
	for (int ii = 0; ii < this->ave_freq.nvars; ii++) {
		char *varname = this->ave_freq.vars[ii].name;
		char *cdfname = this->ave_freq.vars[ii].vname;
		char *varpath = OGS::writeOGSPath(this->ave_freq.vars[ii].path,
			this->timeStepInfo.datetime[ii_tstep],"*");
		// Test if the variable has been activated
		if (this->GetAveFreqArrayStatus(varname)) {
			double *array = NetCDF::readNetCDF(varpath,cdfname,nLon*nLat*nLev);
			vtkFloatArray *vtkarray = createVTKscaf(varname,nLon,nLat,nLev,array);
			this->Mesh->GetCellData()->AddArray(vtkarray);
			vtkarray->Delete(); free(array);
		}
		free(varpath);
		this->UpdateProgress(0.625 + ii*(0.375/ave_freq.nvars));		
	}

	/* SET THE DATETIME FIELD ARRAY

		Set the date and time of the file so that it can be used
		for time stepping purposes. This field contains just information
		and cannot be deactivated.

	*/
	vtkStringArray *strf = createVTKstrf("Date",this->timeStepInfo.datetime[ii_tstep]);
	this->Mesh->GetFieldData()->AddArray(strf);
	strf->Delete();

	// Set output
	output->ShallowCopy(this->Mesh);
	this->UpdateProgress(1.0);

	// Function exit status successful
	return 1;
}

//----------------------------------------------------------------------------
int vtkOGSReader::RequestInformation(vtkInformation* vtkNotUsed(request),
  vtkInformationVector** vtkNotUsed(inputVector), vtkInformationVector* outputVector) {

  	/* READING THE MASTER FILE

		The master file contains the information of the working directory and the
		paths to the mesh and variable files. It also contains a relationship of all
		the existing variables.

		See OGSmesh for further details.
	*/
	OGS::readOGSFile(this->FileName,this->meshfile,
		&this->ave_phys,&this->ave_freq,&this->timeStepInfo);

	/* SCAN THE PHYSICAL VARIABLES

		Detect which physical variables (AVE_PHYS) are present in the master
		file and list them inside the array selection.

	*/
	for(int ii = 0; ii < this->ave_phys.nvars; ii++) {
		// Add the variable to the array selection
		this->AvePhysDataArraySelection->AddArray(this->ave_phys.vars[ii].name);
		this->NumberOfAvePhysComponents += 1;
		if (std::string(this->ave_phys.vars[ii].name) == "Velocity")
			this->NumberOfAvePhysFields += 3;
		else
			this->NumberOfAvePhysFields += 1;
	}

	/* SCAN THE BIOGEOCHEMICAL VARIABLES

		Detect which biogeochemical variables (AVE_FREQ) are present in the master
		file and list them inside the array selection.

	*/
	for(int ii = 0; ii < this->ave_freq.nvars; ii++) {
		// Add the variable to the array selection
		this->AveFreqDataArraySelection->AddArray(this->ave_freq.vars[ii].name);
		this->NumberOfAveFreqComponents += 1;
		this->NumberOfAveFreqFields += 1;		
	}

	/* SET UP THE TIMESTEP

		Time stepping information is contained inside the master file. Here we set
		the time step and the time step range for paraview.

	*/
	vtkInformation* outInfo = outputVector->GetInformationObject(0);

	// Set the time step value
	double *timeSteps;
	timeSteps = (double*)malloc(this->timeStepInfo.ntsteps*sizeof(double));
	for (int ii = 0; ii < this->timeStepInfo.ntsteps; ii++)
		timeSteps[ii] = ii;
    outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_STEPS(), &timeSteps[0], 
    	this->timeStepInfo.ntsteps);
    
	// Set up the time range
    double timeRange[2];
    timeRange[0] = timeSteps[0];
    timeRange[1] = timeSteps[this->timeStepInfo.ntsteps-1];
	outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_RANGE(), timeRange, 2);

	free(timeSteps);

	return 1;
}

//----------------------------------------------------------------------------
void vtkOGSReader::DeleteMesh() {
	this->Mesh->Delete();
	this->Mesh = NULL;
}

//----------------------------------------------------------------------------
void vtkOGSReader::DisableAllAvePhysArrays()
{
  this->AvePhysDataArraySelection->DisableAllArrays();
}

void vtkOGSReader::EnableAllAvePhysArrays()
{
  this->AvePhysDataArraySelection->EnableAllArrays();
}

int vtkOGSReader::GetNumberOfAvePhysArrays()
{
	return this->AvePhysDataArraySelection->GetNumberOfArrays();
}

const char* vtkOGSReader::GetAvePhysArrayName(int index)
{
	if (index >= (int)this->NumberOfAvePhysComponents || index < 0)
	return NULL;
	else
	return this->AvePhysDataArraySelection->GetArrayName(index);
}

int vtkOGSReader::GetAvePhysArrayIndex(const char* name)
{
	return this->AvePhysDataArraySelection->GetArrayIndex(name);
}

int vtkOGSReader::GetAvePhysArrayStatus(const char* name)
{
	return this->AvePhysDataArraySelection->ArrayIsEnabled(name);
}

void vtkOGSReader::SetAvePhysArrayStatus(const char* name, int status)
{
	if (status)
		this->AvePhysDataArraySelection->EnableArray(name);
	else
		this->AvePhysDataArraySelection->DisableArray(name);

	this->Modified();
}

//----------------------------------------------------------------------------
void vtkOGSReader::DisableAllAveFreqArrays()
{
  this->AveFreqDataArraySelection->DisableAllArrays();
}

void vtkOGSReader::EnableAllAveFreqArrays()
{
  this->AveFreqDataArraySelection->EnableAllArrays();
}

int vtkOGSReader::GetNumberOfAveFreqArrays()
{
	return this->AveFreqDataArraySelection->GetNumberOfArrays();
}

const char* vtkOGSReader::GetAveFreqArrayName(int index)
{
	if (index >= (int)this->NumberOfAveFreqComponents || index < 0)
	return NULL;
	else
	return this->AveFreqDataArraySelection->GetArrayName(index);
}

int vtkOGSReader::GetAveFreqArrayIndex(const char* name)
{
	return this->AveFreqDataArraySelection->GetArrayIndex(name);
}

int vtkOGSReader::GetAveFreqArrayStatus(const char* name)
{
	return this->AveFreqDataArraySelection->ArrayIsEnabled(name);
}

void vtkOGSReader::SetAveFreqArrayStatus(const char* name, int status)
{
	if (status)
		this->AveFreqDataArraySelection->EnableArray(name);
	else
		this->AveFreqDataArraySelection->DisableArray(name);

	this->Modified();
}