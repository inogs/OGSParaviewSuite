/*=========================================================================

  Program:   OGSReader
  Module:    vtkOGSReader.cxx

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkOGSReader.h"

#include "vtkCell.h"
#include "vtkPoints.h"
#include "vtkCallbackCommand.h"
#include "vtkCommand.h"
#include "vtkCellData.h"
#include "vtkDataArraySelection.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkRectilinearGrid.h"
#include "vtkFloatArray.h"
#include "vtkDoubleArray.h"
#include "vtkStringArray.h"

#include "vtkObjectFactory.h"

#include <cstdlib>
#include <cstdio>
#include <string>
#include <cstring>

#include <vtksys/SystemTools.hxx>


#ifdef PARAVIEW_USE_MPI
#include "vtkMultiProcessController.h"
vtkCxxSetObjectMacro(vtkOGSReader, Controller, vtkMultiProcessController);
#endif

vtkStandardNewMacro(vtkOGSReader);

//----------------------------------------------------------------------------

/*
	Macro to set the array precision 
*/
#define VTKARRAY vtkDoubleArray

#include "../_utils/field.h"
#include "../_utils/fieldOperations.hpp"
#include "../_utils/OGS.hpp"
#include "../_utils/vtkFields.hpp"
#include "../_utils/netcdfio.hpp"

//----------------------------------------------------------------------------
vtkOGSReader::vtkOGSReader() {

	this->SetNumberOfInputPorts(0);
	this->SetNumberOfOutputPorts(1);

	this->FileName = NULL;

	this->RMeshMask  = 1;
	this->DepthScale = 1000.;
	this->abort      = 0;

	this->Mesh = vtkRectilinearGrid::New();

	this->MaskDataArraySelection    = vtkDataArraySelection::New();
	this->AvePhysDataArraySelection = vtkDataArraySelection::New();
	this->AveFreqDataArraySelection = vtkDataArraySelection::New();
	this->ForcingDataArraySelection = vtkDataArraySelection::New();
	this->GeneralDataArraySelection = vtkDataArraySelection::New();

	#ifdef PARAVIEW_USE_MPI
		this->Controller = NULL;
		this->SetController(vtkMultiProcessController::GetGlobalController());
	#endif
}

//----------------------------------------------------------------------------
vtkOGSReader::~vtkOGSReader() {

	this->FileName = NULL;
	this->SetFileName(NULL);

	this->MaskDataArraySelection->Delete();
	this->AvePhysDataArraySelection->Delete();
	this->AveFreqDataArraySelection->Delete();
	this->ForcingDataArraySelection->Delete();
	this->GeneralDataArraySelection->Delete();

	this->DeleteMesh();

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

	if (this->abort) return 0;

  	// get the info object
	vtkInformation* outInfo = outputVector->GetInformationObject(0);

	// get the output
	vtkRectilinearGrid *output = 
	vtkRectilinearGrid::SafeDownCast(outInfo->Get(vtkRectilinearGrid::DATA_OBJECT()));

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
	if (this->ogsdata.ntsteps() > 1)
		ii_tstep = (int)(requestedTimeValue);

	this->UpdateProgress(0.25);
	
	VTKARRAY *vtkarray;
	int n_vars_loaded = 0;

	/* READING THE PHYSICAL VARIABLES

		Variables inside AVE_PHYS are read here. Velocity is outputed as
		a vector while the others are scalar arrays. User can select which
		variables to load by using a panel.

	*/
	for (int ii = 0; ii < this->ogsdata.var_n(0); ii++) {
		// Test if the variable has been activated
		if (this->GetAvePhysArrayStatus(this->ogsdata.var_name(0,ii))) {
			if (std::string(this->ogsdata.var_name(0,ii)) == "Velocity") {
				field::Field<double>  veloc = NetCDF::readNetCDF2F3(this->ogsdata.var_path(0,ii,ii_tstep).c_str(),
																	"vozocrtx",
																	"vomecrty",
																	"vovecrtz",
																	this->ogsdata.ncells()
																   );							    								  
				if (veloc.isempty()) { vtkErrorMacro("Cannot read NetCDF <Velocity>!Aborting!"); return 0; }

				// We need to project the velocity field from a face centered grid to a cell centered grid
				field::UVW2T(veloc,
							 this->ogsdata.e1(),
							 this->ogsdata.e2(),
							 this->ogsdata.e3(),
							 this->ogsdata.nlon()-1,
							 this->ogsdata.nlat()-1,
							 this->ogsdata.nlev()-1
							);

				vtkarray = VTKFIELD::createVTKfromField<VTKARRAY,double>(this->ogsdata.var_name(0,ii),veloc);
				this->Mesh->GetCellData()->AddArray(vtkarray);
				vtkarray->Delete();
				
				n_vars_loaded++;
			} else {
				field::Field<double>  f = NetCDF::readNetCDF2F(this->ogsdata.var_path(0,ii,ii_tstep).c_str(),
															   this->ogsdata.var_vname(0,ii),
															   this->ogsdata.ncells()
															   );
				if (f.isempty()) { vtkErrorMacro("Cannot read NetCDF <"<<this->ogsdata.var_vname(0,ii)<<">!Aborting!"); return 0; }
				
				vtkarray = VTKFIELD::createVTKfromField<VTKARRAY,double>(this->ogsdata.var_name(0,ii),f);
				this->Mesh->GetCellData()->AddArray(vtkarray);
				vtkarray->Delete();
				
				n_vars_loaded++;
			}
		} else {
			this->Mesh->GetCellData()->RemoveArray(this->ogsdata.var_name(0,ii));
		}
		this->UpdateProgress(0.25 + ii*(0.175/this->ogsdata.var_n(0)));		
	}

	/* READING THE BIOGEOCHEMICAL VARIABLES

		Variables inside AVE_FREQ are read here. User can select which
		variables to load by using a panel.

	*/
	for (int ii = 0; ii < this->ogsdata.var_n(1); ii++) {
		// Test if the variable has been activated
		if (this->GetAveFreqArrayStatus(this->ogsdata.var_name(1,ii))) {
				field::Field<double>  f = NetCDF::readNetCDF2F(this->ogsdata.var_path(1,ii,ii_tstep).c_str(),
															   this->ogsdata.var_vname(1,ii),
															   this->ogsdata.ncells()
															   );
				if (f.isempty()) { vtkErrorMacro("Cannot read NetCDF <"<<this->ogsdata.var_vname(1,ii)<<">!Aborting!"); return 0; }
				
				vtkarray = VTKFIELD::createVTKfromField<VTKARRAY,double>(this->ogsdata.var_name(1,ii),f);
				this->Mesh->GetCellData()->AddArray(vtkarray);
				vtkarray->Delete();
				
				n_vars_loaded++;
		} else {
			this->Mesh->GetCellData()->RemoveArray(this->ogsdata.var_name(1,ii));
		}
		this->UpdateProgress(0.425 + ii*(0.175/this->ogsdata.var_n(1)));		
	}

	/* READING THE FORCINGS VARIABLES

		Variables inside FORCINCS are read here. User can select which
		variables to load by using a panel.
	*/
	for (int ii = 0; ii < this->ogsdata.var_n(2); ii++) {
		// Test if the variable has been activated
		if (this->GetForcingArrayStatus(this->ogsdata.var_name(2,ii))) {
				field::Field<double>  f = NetCDF::readNetCDF2F(this->ogsdata.var_path(2,ii,ii_tstep).c_str(),
															   this->ogsdata.var_vname(2,ii),
															   this->ogsdata.ncells()
															   );
				if (f.isempty()) { vtkErrorMacro("Cannot read NetCDF <"<<this->ogsdata.var_vname(2,ii)<<">!Aborting!"); return 0; }
				
				vtkarray = VTKFIELD::createVTKfromField<VTKARRAY,double>(this->ogsdata.var_name(2,ii),f);
				this->Mesh->GetCellData()->AddArray(vtkarray);
				vtkarray->Delete();
				
				n_vars_loaded++;
		} else {
			this->Mesh->GetCellData()->RemoveArray(this->ogsdata.var_name(2,ii));
		}
		this->UpdateProgress(0.6 + ii*(0.175/this->ogsdata.var_n(2)));		
	}

	/* READING THE GENERAL VARIABLES
		Variables inside GENERAL are read here. User can select which
		variables to load by using a panel.
	*/
	for (int ii = 0; ii < this->ogsdata.var_n(3); ii++) {
		// Test if the variable has been activated
		if (this->GetGeneralArrayStatus(this->ogsdata.var_name(3,ii))) {
				field::Field<double>  f = NetCDF::readNetCDF2F(this->ogsdata.var_path(3,ii,ii_tstep).c_str(),
															   this->ogsdata.var_vname(3,ii),
															   this->ogsdata.ncells()
															   );
				if (f.isempty()) { vtkErrorMacro("Cannot read NetCDF <"<<this->ogsdata.var_vname(3,ii)<<">!Aborting!"); return 0; }
				
				vtkarray = VTKFIELD::createVTKfromField<VTKARRAY,double>(this->ogsdata.var_name(3,ii),f);
				this->Mesh->GetCellData()->AddArray(vtkarray);
				vtkarray->Delete();
				
				n_vars_loaded++;
		} else {
			this->Mesh->GetCellData()->RemoveArray(this->ogsdata.var_name(3,ii));
		}
		this->UpdateProgress(0.775 + ii*(0.175/this->ogsdata.var_n(3)));		
	}

	/* SET THE METADATA ARRAY
		Set the metadata, an array that contains multiple information for
		the performance of the OGS filters. Data:
			0 -> Date
			1 -> Datevec
			2 -> Conversion factors
			3 -> Loaded variables
			4 -> File name
			5 -> Meshmask
	*/
	std::string aux_str;
	vtkStringArray *vtkmetadata = VTKFIELD::createVTKstrf("Metadata",6,NULL);

	// Set the current file date
	vtkmetadata->SetValue(0,this->ogsdata.datetime(ii_tstep));

	// Set the datevec
	aux_str = std::to_string(this->ogsdata.ntsteps()) + std::string(";");
	for (int ii = 0; ii < this->ogsdata.ntsteps(); ii++)
		aux_str += std::string(this->ogsdata.datetime(ii)) + std::string(";");
	vtkmetadata->SetValue(1,aux_str.c_str());

	// Set conversion factors
	aux_str = std::to_string(this->DepthScale);
	vtkmetadata->SetValue(2,aux_str.c_str());

	// Set the number of variables loaded
	aux_str = std::to_string(n_vars_loaded) + std::string(";");
	// AVE_PHYS
	for (int ii = 0; ii < this->ogsdata.var_n(0); ii++) 
		if (this->GetAvePhysArrayStatus(this->ogsdata.var_name(0,ii)))
			aux_str += std::string(this->ogsdata.var_name(0,ii)) + std::string(";");
	// AVE_FREQ
	for (int ii = 0; ii < this->ogsdata.var_n(1); ii++) 
		if (this->GetAveFreqArrayStatus(this->ogsdata.var_name(1,ii)))
			aux_str += std::string(this->ogsdata.var_name(1,ii)) + std::string(";");
	// FORCINGS
	for (int ii = 0; ii < this->ogsdata.var_n(2); ii++) 
		if (this->GetForcingArrayStatus(this->ogsdata.var_name(2,ii)))
			aux_str += std::string(this->ogsdata.var_name(2,ii)) + std::string(";");
	// GENERALS
	for (int ii = 0; ii < this->ogsdata.var_n(3); ii++) 
		if (this->GetGeneralArrayStatus(this->ogsdata.var_name(3,ii)))
			aux_str += std::string(this->ogsdata.var_name(3,ii)) + std::string(";");
	vtkmetadata->SetValue(3,aux_str.c_str());

	// Set the file name
	vtkmetadata->SetValue(4,this->FileName);

	// Set the mesh name
	aux_str = this->ogsdata.meshmask() + std::string(";") + this->ogsdata.meshfile();
	vtkmetadata->SetValue(5,aux_str.c_str());

	// Add array to mesh
	this->Mesh->GetFieldData()->AddArray(vtkmetadata);
	vtkmetadata->Delete();

	/* 
		FINAL TOUCHES
	*/
	output->ShallowCopy(this->Mesh);
	this->UpdateProgress(1.0);

	// Function exit status successful
	return 1;
}

//----------------------------------------------------------------------------
int vtkOGSReader::RequestInformation(vtkInformation* vtkNotUsed(request),
  vtkInformationVector** vtkNotUsed(inputVector), vtkInformationVector* outputVector) {

  	vtkInformation* outInfo = outputVector->GetInformationObject(0);

  	/* SET UP THE MASKS

  		 Set up the data array selection containing the masks to be loaded.

  	*/
  	this->MaskDataArraySelection->AddArray("Sub-basins");         // Mask containing the different sub-basins of the MED
  	this->MaskDataArraySelection->AddArray("Continental shelf");  // Mask of the continental shelf

  	this->UpdateProgress(0.0);

  	/* READING THE MASTER FILE

		The master file contains the information of the working directory and the
		paths to the mesh and variable files. It also contains a relationship of all
		the existing variables.

		See OGSmesh for further details.
	*/
	this->ogsdata.SetFile(this->FileName);
	if (ogsdata.readMainFile() == -1) { 
		vtkErrorMacro("Cannot read <"<<this->FileName<<">!\nAborting");
		return 0; this->abort = 1;
	}

	/* READING THE MESH FILE 
		
		The mesh data (dims, Lon2Meters, Lat2Meters and nav_lev) are inside
		the .ogsmsh file containing the correct resolution of the meshmask for
		the current simulation.

		This is done in the RequestInformation to ensure it is only read once
		and stored in memory.

		See OGS.hpp/OGS.cpp for further details.

	*/
	if (this->ogsdata.readMesh() < 0) {
		vtkErrorMacro("Problems reading the mesh!\nAborting.");
		return 0; this->abort = 1;		
	}

	this->UpdateProgress(0.05);

	/* CREATE RECTILINEAR GRID

		The rectilinear grid is created here from the mesh data that has
		been previously read. This is handled in the RequestInformation as it only needs to
		be executed once. Moreover, the rectilinear grid is only created when the mesh 
		dimensions are (0,0,0).

	*/
	VTKFIELD::createRectilinearGrid(this->ogsdata.nlon(),
									this->ogsdata.nlat(),
									this->ogsdata.nlev(),
									this->ogsdata.lon2meters(),
									this->ogsdata.lat2meters(),
									this->ogsdata.nav_lev(),
									this->DepthScale,
									this->Mesh
								   );

	this->UpdateProgress(0.10);

	/* ADD MASK ARRAYS

		Add the masks to the mesh. This process is handled here since it only needs to change when
		the request information executes, and not with the time-steppin

		Current masks are:
			> Sub-basins: mask with all the basins of the MED.
			> Continental shelf: area from 0 to 200 m of depth.

	*/
	VTKARRAY *vtkarray;

	if (this->GetMaskArrayStatus("Sub-basins")) {
		vtkarray = VTKFIELD::createVTKfromField<VTKARRAY,double>("basins mask",this->ogsdata.mask(0));
		this->Mesh->GetCellData()->AddArray(vtkarray);
		vtkarray->Delete();
	} else {
		this->Mesh->GetCellData()->RemoveArray("basins mask");
	}

	// Continental shelf mask ("coast mask")
	if (this->GetMaskArrayStatus("Continental shelf")) {
		vtkarray = VTKFIELD::createVTKfromField<VTKARRAY,double>("coast mask",this->ogsdata.mask(1));
		this->Mesh->GetCellData()->AddArray(vtkarray);
		vtkarray->Delete();
	} else {
		this->Mesh->GetCellData()->RemoveArray("coast mask");
	}

	this->UpdateProgress(0.15);

	/* ADD MESHMASK STRETCHING ARRAYS

		Add the stretching arrays e1, e2 and e3 found in the meshmask. These arrays are needed to
		project the velocity from a face centered coordinates to cell centered coordinates as well
		as to compute gradients.

	*/
	this->ogsdata.readMeshmask();

	if (this->RMeshMask) {
		// e1
		vtkarray = VTKFIELD::createVTKfromField<VTKARRAY,double>("e1",this->ogsdata.e1());
		this->Mesh->GetCellData()->AddArray(vtkarray);
		vtkarray->Delete();
		// e2
		vtkarray = VTKFIELD::createVTKfromField<VTKARRAY,double>("e2",this->ogsdata.e2());
		this->Mesh->GetCellData()->AddArray(vtkarray);
		vtkarray->Delete();
		// e3
		vtkarray = VTKFIELD::createVTKfromField<VTKARRAY,double>("e3",this->ogsdata.e3());
		this->Mesh->GetCellData()->AddArray(vtkarray);
		vtkarray->Delete();
	} else {
		this->Mesh->GetCellData()->RemoveArray("e1");
		this->Mesh->GetCellData()->RemoveArray("e2");
		this->Mesh->GetCellData()->RemoveArray("e3");
	}

	this->UpdateProgress(0.20);

	/* SCAN THE PHYSICAL VARIABLES

		Detect which physical variables (AVE_PHYS) are present in the master
		file and list them inside the array selection.

	*/
	for(int ii = 0; ii < this->ogsdata.var_n(0); ii++)
		this->AvePhysDataArraySelection->AddArray(this->ogsdata.var_name(0,ii));

	/* SCAN THE BIOGEOCHEMICAL VARIABLES

		Detect which biogeochemical variables (AVE_FREQ) are present in the master
		file and list them inside the array selection.

	*/
	for(int ii = 0; ii < this->ogsdata.var_n(1); ii++)
		this->AveFreqDataArraySelection->AddArray(this->ogsdata.var_name(1,ii));

	/* SCAN THE FORCING VARIABLES

		Detect which forcing variables (FORCINGS) are present in the master
		file and list them inside the array selection.

	*/
	for(int ii = 0; ii < this->ogsdata.var_n(2); ii++)
		this->ForcingDataArraySelection->AddArray(this->ogsdata.var_name(2,ii));

	/* SCAN THE GENERAL VARIABLES

		Detect which general variables (GENERAL) are present in the master
		file and list them inside the array selection.

	*/
	for(int ii = 0; ii < this->ogsdata.var_n(3); ii++)
		this->GeneralDataArraySelection->AddArray(this->ogsdata.var_name(3,ii));


	/* SET UP THE TIMESTEP

		Time stepping information is contained inside the master file. Here we set
		the time step and the time step range for paraview.

	*/
	// Set the time step value
	double *timeSteps = NULL; timeSteps = new double[this->ogsdata.ntsteps()];
	for (int ii = 0; ii < this->ogsdata.ntsteps(); ii++)
		timeSteps[ii] = ii;
    outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_STEPS(), 
    	&timeSteps[0], this->ogsdata.ntsteps());
    
	// Set up the time range
    double timeRange[2] = {timeSteps[0], timeSteps[this->ogsdata.ntsteps()-1]};
	outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_RANGE(), timeRange, 2);

	delete [] timeSteps;

	this->UpdateProgress(0.25);
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
	if (index >= (int)this->GetNumberOfAvePhysArrays() || index < 0)
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
	if (index >= (int)this->GetNumberOfAveFreqArrays() || index < 0)
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

//----------------------------------------------------------------------------
void vtkOGSReader::DisableAllMaskArrays()
{
  this->MaskDataArraySelection->DisableAllArrays();
}

void vtkOGSReader::EnableAllMaskArrays()
{
  this->MaskDataArraySelection->EnableAllArrays();
}

int vtkOGSReader::GetNumberOfMaskArrays()
{
	return this->MaskDataArraySelection->GetNumberOfArrays();
}

const char* vtkOGSReader::GetMaskArrayName(int index)
{
	if (index >= (int)this->GetNumberOfMaskArrays() || index < 0)
		return NULL;
	else
		return this->MaskDataArraySelection->GetArrayName(index);
}

int vtkOGSReader::GetMaskArrayIndex(const char* name)
{
	return this->MaskDataArraySelection->GetArrayIndex(name);
}

int vtkOGSReader::GetMaskArrayStatus(const char* name)
{
	return this->MaskDataArraySelection->ArrayIsEnabled(name);
}

void vtkOGSReader::SetMaskArrayStatus(const char* name, int status)
{
	if (status)
		this->MaskDataArraySelection->EnableArray(name);
	else
		this->MaskDataArraySelection->DisableArray(name);

	this->Modified();
}

//----------------------------------------------------------------------------
void vtkOGSReader::DisableAllForcingArrays()
{
  this->ForcingDataArraySelection->DisableAllArrays();
}

void vtkOGSReader::EnableAllForcingArrays()
{
  this->ForcingDataArraySelection->EnableAllArrays();
}

int vtkOGSReader::GetNumberOfForcingArrays()
{
	return this->ForcingDataArraySelection->GetNumberOfArrays();
}

const char* vtkOGSReader::GetForcingArrayName(int index)
{
	if (index >= (int)this->GetNumberOfForcingArrays() || index < 0)
		return NULL;
	else
		return this->ForcingDataArraySelection->GetArrayName(index);
}

int vtkOGSReader::GetForcingArrayIndex(const char* name)
{
	return this->ForcingDataArraySelection->GetArrayIndex(name);
}

int vtkOGSReader::GetForcingArrayStatus(const char* name)
{
	return this->ForcingDataArraySelection->ArrayIsEnabled(name);
}

void vtkOGSReader::SetForcingArrayStatus(const char* name, int status)
{
	if (status)
		this->ForcingDataArraySelection->EnableArray(name);
	else
		this->ForcingDataArraySelection->DisableArray(name);

	this->Modified();
}

//----------------------------------------------------------------------------
void vtkOGSReader::DisableAllGeneralArrays()
{
  this->GeneralDataArraySelection->DisableAllArrays();
}

void vtkOGSReader::EnableAllGeneralArrays()
{
  this->GeneralDataArraySelection->EnableAllArrays();
}

int vtkOGSReader::GetNumberOfGeneralArrays()
{
	return this->GeneralDataArraySelection->GetNumberOfArrays();
}

const char* vtkOGSReader::GetGeneralArrayName(int index)
{
	if (index >= (int)this->GetNumberOfGeneralArrays() || index < 0)
		return NULL;
	else
		return this->GeneralDataArraySelection->GetArrayName(index);
}

int vtkOGSReader::GetGeneralArrayIndex(const char* name)
{
	return this->GeneralDataArraySelection->GetArrayIndex(name);
}

int vtkOGSReader::GetGeneralArrayStatus(const char* name)
{
	return this->GeneralDataArraySelection->ArrayIsEnabled(name);
}

void vtkOGSReader::SetGeneralArrayStatus(const char* name, int status)
{
	if (status)
		this->GeneralDataArraySelection->EnableArray(name);
	else
		this->GeneralDataArraySelection->DisableArray(name);

	this->Modified();
}