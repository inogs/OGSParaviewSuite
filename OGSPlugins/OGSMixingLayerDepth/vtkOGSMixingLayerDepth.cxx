/*=========================================================================

  Program:   OGSMixingLayerDepth
  Module:    vtkOGSMixingLayerDepth.cxx

  Copyright (c) 2020 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkIntArray.h"
#include "vtkFloatArray.h"
#include "vtkDoubleArray.h"
#include "vtkTypeUInt8Array.h"
#include "vtkStringArray.h"
#include "vtkCell.h"
#include "vtkPoints.h"
#include "vtkCellData.h"
#include "vtkPointData.h"
#include "vtkFieldData.h"
#include "vtkDataSet.h"
#include "vtkRectilinearGrid.h"
#include "vtkDataArraySelection.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#include "vtkOGSMixingLayerDepth.h"

#include "vtkObjectFactory.h"

#include <string>
#include <map>
#include <algorithm>

#ifdef PARAVIEW_USE_MPI
#include "vtkMultiProcessController.h"
vtkCxxSetObjectMacro(vtkOGSMixingLayerDepth, Controller, vtkMultiProcessController);
#endif

vtkStandardNewMacro(vtkOGSMixingLayerDepth);

//----------------------------------------------------------------------------
#include "macros.h"
#include "fieldOperations.h"
#include "vtkFields.h"
#include "vtkOperations.h"

//----------------------------------------------------------------------------
vtkOGSMixingLayerDepth::vtkOGSMixingLayerDepth() {
	this->useDensity = false;
	this->zref       = 10.;
	this->dT         = 0.2;
	this->drho       = 0.03;
	this->mask_field = NULL;
	this->epsi       = 0.001;
	this->iscelld    = true;
	this->isReqInfo  = false;
	this->nProcs     = 0;
	this->procId     = 0;

	#ifdef PARAVIEW_USE_MPI
		this->Controller = NULL;
		this->SetController(vtkMultiProcessController::GetGlobalController());
	#endif
}

//----------------------------------------------------------------------------
vtkOGSMixingLayerDepth::~vtkOGSMixingLayerDepth() {
	this->Setmask_field(NULL);
	
	#ifdef PARAVIEW_USE_MPI
		this->SetController(NULL);	
	#endif
}

//----------------------------------------------------------------------------
int vtkOGSMixingLayerDepth::RequestInformation(vtkInformation* vtkNotUsed(request),
  vtkInformationVector** vtkNotUsed(inputVector), vtkInformationVector* outputVector) {

	/* SET UP THE PARALLEL CONTROLLER

		The MPI threads come initialized by the ParaView server. Here
		we set up the environment for this filter.

	*/
	#ifdef PARAVIEW_USE_MPI
	if (this->Controller->GetNumberOfProcesses() > 1) {
		this->nProcs = this->Controller->GetNumberOfProcesses();
		this->procId = this->Controller->GetLocalProcessId();
	}

	// Stop all threads except from the master to execute
	if (this->procId > 0) return 1;
	#endif

  	this->isReqInfo = true;
  	return 1;
}

//----------------------------------------------------------------------------
int vtkOGSMixingLayerDepth::RequestData(vtkInformation *vtkNotUsed(request),
	vtkInformationVector **inputVector,vtkInformationVector *outputVector) {
	// Get the info objects
	vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
	vtkInformation *outInfo = outputVector->GetInformationObject(0);

	// Stop all threads except from the master to execute
	#ifdef PARAVIEW_USE_MPI
	if (this->procId > 0) return 1;
	#endif

	// Get the input and output
	vtkRectilinearGrid *input = vtkRectilinearGrid::SafeDownCast(
		inInfo->Get(vtkDataObject::DATA_OBJECT()));
	vtkRectilinearGrid *output = vtkRectilinearGrid::SafeDownCast(
		outInfo->Get(vtkDataObject::DATA_OBJECT()));

	// We just want to copy the mesh, not the variables
	output->ShallowCopy(input);
	output->GetFieldData()->PassData(input->GetFieldData());

	this->UpdateProgress(0.);

	// Test for cell data or point data
	if (VTKARRAY::SafeDownCast(input->GetCellData()->GetArray("e1")) == NULL) {
		this->iscelld  = false; 
	}

	// This section is only executed once, to populate the xyz and
	// cId2zId arrays. Successive iterations should not execute.
	// This section is included here since RequestInformation gives
	// troubles when restarting.
	if (this->xyz.isempty() || this->isReqInfo) {

		this->isReqInfo = false;
		
		// Recover Metadata array (depth factor)
		vtkStringArray *vtkmetadata = vtkStringArray::SafeDownCast(
			input->GetFieldData()->GetAbstractArray("Metadata"));
		double dfact = (vtkmetadata != NULL) ? std::stod( vtkmetadata->GetValue(2) ) : 1000.;
		
		if (vtkmetadata == NULL) 
			vtkWarningMacro("Field array Metadata not found! Depth factor set to 1000. automatically.");

		// Recover cell or point coordinates, depending on the array type
		this->xyz = (this->iscelld) ? VTK::getVTKCellCenters(input,dfact) : 
		                              VTK::getVTKCellPoints(input,dfact);

		// Up to this point we have the cell centers or point coordinates correctly
		// stored under "xyz". Now we shall find the number of unique z coordinates or,
		// depending on the user input, the coordinates of each depth level, as well as
		// its mesh connectivity (cId2zId).
		this->zcoords.clear();
		this->cId2zId = field::countDepthLevels(this->xyz,this->zcoords,this->epsi,true);
	}

	this->UpdateProgress(.1);

	// Output connectivity (Uncomment for debugging purposes)
	//vtkIntArray *vtkcId2zId;
	//vtkcId2zId = VTK::createVTKfromField<vtkIntArray,int>("cId2zId",this->cId2zId);
	//output->GetCellData()->AddArray(vtkcId2zId);
	//vtkcId2zId->Delete();

	// Coast mask
	VTKMASK *vtkMask;
	vtkMask = VTKMASK::SafeDownCast( 
		(this->iscelld) ? input->GetCellData()->GetArray(this->mask_field) : 
						  input->GetPointData()->GetArray(this->mask_field) );
	if (!vtkMask) {
		vtkErrorMacro("Error loading <"<<this->mask_field<<">! Aborting!");
		return 0;
	}
	field::Field<FLDMASK> cmask;
	cmask = VTK::createFieldfromVTK<VTKMASK,FLDMASK>(vtkMask);	

	// Compute MLD according to the definition
	field::Field<FLDARRAY> mld(this->xyz.len(),1,0.);
	field::Field<FLDMASK>  mldmask(this->xyz.len(),1,(FLDMASK)(0));

	// Load array, either density or temperature
	VTKARRAY *vtkArray; double dMLD;
	if (this->useDensity) {
		dMLD = this->drho;
		// Recover array
		vtkArray = VTKARRAY::SafeDownCast( 
			(this->iscelld) ? input->GetCellData()->GetArray("Density") : 
							  input->GetPointData()->GetArray("Density") );	
	} else {
		dMLD = this->dT;
		// Recover array
		vtkArray = VTKARRAY::SafeDownCast( 
			(this->iscelld) ? input->GetCellData()->GetArray("Temperature") : 
							  input->GetPointData()->GetArray("Temperature") );	
	}
	if (!vtkArray) {
		vtkErrorMacro("Error loading Temperature or Density! Aborting!");
		return 0;
	}
	field::Field<FLDARRAY> array;
	array = VTK::createFieldfromVTK<VTKARRAY,FLDARRAY>(vtkArray);

	// First find the reference depth value
	double diff = 1e99, diffold = 1e99; int zIdRef = -1;
	for (int zId = 0; zId < this->zcoords.size(); ++zId) {
		diff = (double)(fabs(this->zref + zcoords[zId])); // zcoords is negative!!
		if (diff < diffold) { zIdRef = zId; diffold = diff; }
	}
	if (zIdRef < 0) { 
		vtkErrorMacro("Problems finding reference depth! Aborting!");
		return 0;
	}
	this->UpdateProgress(0.25);

	// Create a map to refer the data on the layers to the data on the mesh.
	// Populate them by looping the mesh once.
	std::map<int,std::vector<int>> mMeshPerLayer;

	// Set up the different layers
	#pragma omp parallel
	{
	// Set up the different layers
	for (int ii=OMP_THREAD_NUM; ii < array.get_n(); ii+=OMP_NUM_THREADS) {
		// Set maps
		int ind  = cId2zId[ii][0];
		#pragma omp critical
		{
		mMeshPerLayer[ind].push_back(ii);
		}	
	}
	}
	this->UpdateProgress(0.5);

	// For each depth layer, compute the mask and store MLD
	// This algorithm requires to work on a rectilinear grid with
	// the same number of points at each plane
	for (int zId = 0; zId < this->zcoords.size(); zId += 1) {
		// Discard all that is above the reference
		if (zId <= zIdRef) {
			for (int ii=0; ii<mMeshPerLayer[zId].size(); ++ii) {
				int ind    = mMeshPerLayer[zId][ii];
				int ind0   = mMeshPerLayer[0][ii];

				if (cmask[ind][0] == 0) continue;

				mld[ind0][0]    = -zcoords[zIdRef]; // zcoords is negative!!
				mldmask[ind][0] = 1;
			}
			continue;
		}
		for (int ii=0; ii<mMeshPerLayer[zId].size(); ++ii) {
			// Indices
			int ind    = mMeshPerLayer[zId][ii];
			int ind0   = mMeshPerLayer[0][ii];
			int ind1   = mMeshPerLayer[zId-1][ii];
			int indref = mMeshPerLayer[zIdRef][ii];

			double Val  = array[ind][0];
			double Vref = array[indref][0];

			if (cmask[ind][0] == 0) continue;

			// Compute MLD
			if (mldmask[ind1][0] && (double)(fabs(Val-Vref)) < dMLD) {
				mld[ind0][0]    = -zcoords[zId]; // zcoords is negative!!
				mldmask[ind][0] = 1;	
			}
		}
	}
	this->UpdateProgress(0.65);
	
	// Expand MLD to all the layers
	for (int zId = 1; zId < this->zcoords.size(); zId += 1) {
		for (int ii=0; ii<mMeshPerLayer[zId].size(); ++ii) {
			int ind0   = mMeshPerLayer[0][ii];
			int ind    = mMeshPerLayer[zId][ii];

			if (cmask[ind][0] == 0) continue;

			mld[ind][0] = mld[ind0][0];
		}
	}
	this->UpdateProgress(0.75);

	// Output field
	VTKARRAY *vtkmld; VTKMASK *vtkmldmask;
	vtkmld     = VTK::createVTKfromField<VTKARRAY,FLDARRAY>("Mixing Layer Depth",mld);
	vtkmldmask = VTK::createVTKfromField<VTKMASK,FLDMASK>("Mixing Layer Depth Mask",mldmask);

	if (iscelld) {
		output->GetCellData()->AddArray(vtkmld);
		output->GetCellData()->AddArray(vtkmldmask);
	} else {
		output->GetPointData()->AddArray(vtkmld);
		output->GetPointData()->AddArray(vtkmldmask);
	}

	vtkmld->Delete();
	vtkmldmask->Delete();

	// Return
	this->UpdateProgress(1.0);
	return 1;
}
