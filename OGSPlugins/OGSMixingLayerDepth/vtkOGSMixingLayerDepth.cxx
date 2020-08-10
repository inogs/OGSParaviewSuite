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
	std::vector<int> mMeshSurface;
	std::map<int,std::vector<int>> mMeshPerDepth;

	// Set up a connectivity matrix that relates each point of the surface with
	// the points that have same x, y coordinates but different z.
	// Loop the mesh
	for (int ii=0; ii < this->xyz.len(); ii+=1) {
		// is the point on the coast? if so discard it
		if (cmask[ii][0] == 0) continue;
		v3::V3 p1 = this->xyz[ii];
		double p1z = fabs(p1[2]); p1[2] = 0.; // So that we can do the norm later
		// is the point on the surface? if not discard it
		if (fabs(p1z - fabs(this->zcoords[0])) > 1e-3) continue;

		// At this point we should have a valid point and we know
		// how many depth levels there are. We will loop on the depth levels
		// and find the cell that is right under the current one.
		mMeshSurface.push_back(ii); // Points that are on the surface
		mMeshPerDepth[ii].push_back(ii);
		int cId = ii; 
		for (int zId = 1; zId < this->zcoords.size(); zId += 1) {
			// Find the cell neighbours
			std::vector<int> neighbours = (this->iscelld) ? 
				VTK::getConnectedCells(output,cId) : VTK::getConnectedPoints(output,cId);
			// From all the neighbours find the one that has the
			// same x, y and greater z
			int nId = -1; double diff = 1e99;
			for (int jj = 0; jj < neighbours.size(); ++jj) {
				int ind = neighbours[jj];
				v3::V3 p2 = this->xyz[ind];
				double p2z = fabs(p2[2]); p2[2] = 0.;

				if (cmask[ind][0] == 0) continue;
				if (fabs(p2z - fabs(this->zcoords[zId])) > 1e-3) continue;

				double d = (p1 - p2).norm2();
				if (d < diff) { diff = d; nId = jj; }
			}
			if (nId < 0) break; // No more depths have been found
			// Store connectivity
			mMeshPerDepth[ii].push_back(neighbours[nId]); 
			cId = neighbours[nId]; // Find next cell
		}
	}
	this->UpdateProgress(0.5);

	// Loop the points that are on the surface
	for (int ii=0; ii < mMeshSurface.size(); ++ii) {
		int ind   = mMeshSurface[ii];
		int ziref = mMeshPerDepth[ind][zIdRef]; // Reference indices
		// Loop the points that are on the water column
		// Find argmin(fabs(array[z]-array[zref])) < darray
		int idZmin = -1;
		for (int jj = 0; jj<mMeshPerDepth[ind].size(); ++jj) {
			int zind = mMeshPerDepth[ind][jj];
			
			double Val  = array[zind][0];
			double Vref = array[ziref][0];

			// Compute MLD
			if ((double)(fabs(Val-Vref)) < dMLD) idZmin = jj;
		}
		if (idZmin < 0) idZmin = zIdRef;
		// Write MLD and MLD mask arrays
		for (int jj = 0; jj<mMeshPerDepth[ind].size(); ++jj) {
			int zind = mMeshPerDepth[ind][jj];

			mld[zind][0] = -zcoords[idZmin]; // zcoords is negative!!
			if (jj < idZmin) mldmask[zind][0] = 1;
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
