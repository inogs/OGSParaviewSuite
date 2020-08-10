/*=========================================================================

  Program:   OGSRossbyRadius
  Module:    vtkOGSRossbyRadius.cxx

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

#include "vtkOGSRossbyRadius.h"

#include "vtkObjectFactory.h"

#include <string>
#include <map>
#include <algorithm>

#ifdef PARAVIEW_USE_MPI
#include "vtkMultiProcessController.h"
vtkCxxSetObjectMacro(vtkOGSRossbyRadius, Controller, vtkMultiProcessController);
#endif

vtkStandardNewMacro(vtkOGSRossbyRadius);

//----------------------------------------------------------------------------
#include "macros.h"
#include "projection.h"
#include "fieldOperations.h"
#include "vtkFields.h"
#include "vtkOperations.h"

//----------------------------------------------------------------------------
vtkOGSRossbyRadius::vtkOGSRossbyRadius() {
	this->g          = 9.81;
	this->f_cor_ct   = 1.0e-4;
	this->Omega      = 7.2921e-5;
	this->useCtFcor  = false;

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
vtkOGSRossbyRadius::~vtkOGSRossbyRadius() {
	this->Setmask_field(NULL);
	
	#ifdef PARAVIEW_USE_MPI
		this->SetController(NULL);	
	#endif
}

//----------------------------------------------------------------------------
int vtkOGSRossbyRadius::RequestInformation(vtkInformation* vtkNotUsed(request),
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
int vtkOGSRossbyRadius::RequestData(vtkInformation *vtkNotUsed(request),
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
	std::string projName;
	vtkStringArray *vtkmetadata = vtkStringArray::SafeDownCast(
		input->GetFieldData()->GetAbstractArray("Metadata"));
	
	projName = (vtkmetadata != NULL) ? vtkmetadata->GetValue(7) : std::string("Mercator");
	std::transform(projName.begin(), projName.end(), projName.begin(), ::tolower);

	if (this->xyz.isempty() || this->isReqInfo) {

		this->isReqInfo = false;
		
		// Recover Metadata array (depth factor)

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

	// e3
	VTKARRAY *vtke3;
	vtke3 = VTKARRAY::SafeDownCast( 
		(this->iscelld) ? input->GetCellData()->GetArray("e3") : 
						  input->GetPointData()->GetArray("e3") );
	if (!vtke3) {
		vtkErrorMacro("Error loading <e3>! Aborting!");
		return 0;
	}
	field::Field<FLDARRAY> e3;
	e3 = VTK::createFieldfromVTK<VTKARRAY,FLDARRAY>(vtke3);

	// Density
	VTKARRAY *vtkrho;
	vtkrho = VTKARRAY::SafeDownCast( 
		(this->iscelld) ? input->GetCellData()->GetArray("Density") : 
						  input->GetPointData()->GetArray("Density") );
	if (!vtkrho) {
		vtkErrorMacro("Error loading <Density>! Aborting!");
		return 0;
	}	
	field::Field<FLDARRAY> rho;
	rho = VTK::createFieldfromVTK<VTKARRAY,FLDARRAY>(vtkrho);

	// MLDmask mask
	VTKMASK *vtkmldmask;
	vtkmldmask = VTKMASK::SafeDownCast( 
		(this->iscelld) ? input->GetCellData()->GetArray("Mixing Layer Depth Mask") : 
						  input->GetPointData()->GetArray("Mixing Layer Depth Mask") );
	if (!vtkmldmask) {
		vtkErrorMacro("Error loading <Mixing Layer Depth Mask>! Aborting!");
		return 0;
	}
	field::Field<FLDMASK> mldmask;
	mldmask = VTK::createFieldfromVTK<VTKMASK,FLDMASK>(vtkmldmask);

	// MLD
	VTKARRAY *vtkmld;
	vtkmld = VTKARRAY::SafeDownCast( 
		(this->iscelld) ? input->GetCellData()->GetArray("Mixing Layer Depth") : 
						  input->GetPointData()->GetArray("Mixing Layer Depth") );
	if (!vtkmld) {
		vtkErrorMacro("Error loading <Mixing Layer Depth>! Aborting!");
		return 0;
	}	
	field::Field<FLDARRAY> mld;
	mld = VTK::createFieldfromVTK<VTKARRAY,FLDARRAY>(vtkmld);

	// Computation arrays
	field::Field<FLDARRAY> r_ext(this->xyz.len(),1,0.);
	field::Field<FLDARRAY> r_int(this->xyz.len(),1,0.);
	field::Field<FLDARRAY> rho1_aux(this->xyz.len(),1,0.);
	field::Field<FLDARRAY> rho2_aux(this->xyz.len(),1,0.);
	field::Field<FLDARRAY> d1_aux(this->xyz.len(),1,0.);
	field::Field<FLDARRAY> d2_aux(this->xyz.len(),1,0.);
	field::Field<FLDARRAY> f_cor(this->xyz.len(),1,0.);

	this->UpdateProgress(0.25);

	// Create a map to refer the data on the layers to the data on the mesh.
	// Populate them by looping the mesh once.
	std::vector<int> mMeshSurface;
	std::map<int,std::vector<int>> mMeshPerDepth;

	// Projection
	PROJ::Projection p;

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

				//if (cmask[ind][0] == 0) continue;
				if (fabs(p2z - fabs(this->zcoords[zId])) > 1e-3) continue;

				double d = (p1 - p2).norm2();
				if (d < diff) { diff = d; nId = jj; }
			}
			if (nId < 0) break; // No more depths have been found
			// Store connectivity
			mMeshPerDepth[ii].push_back(neighbours[nId]); 
			cId = neighbours[nId]; // Find next cell
		}

		// Obtain f_cor
		if (this->useCtFcor) {
			f_cor[ii][0] = this->f_cor_ct;
		} else {
			double lon = this->xyz[ii][0], lat = this->xyz[ii][1];
			p.transform_point(projName,"degrees",lon,lat);
			f_cor[ii][0] = 2*this->Omega*sin(PROJ::deg2rad(lat));
		}
	}
	this->UpdateProgress(0.5);

	// Loop the points that are on the surface
	for (int ii=0; ii < mMeshSurface.size(); ++ii) {
		int ind = mMeshSurface[ii];

		// First loop at the depth levels
		// Compute Barotropic Rossby radius and density averages
		for (int jj = 0; jj<mMeshPerDepth[ind].size(); ++jj) {
			int zind = mMeshPerDepth[ind][jj];

			if (cmask[zind][0] == 0) continue;

			// Barotropic Rossby radius of deformation (external)
			// LR = sqrt(g*D)/f_cor
			r_ext[ind][0] = sqrt(this->g*-this->zcoords[jj])/f_cor[ind][0]/1e3; // [km] zcoords is negative
			
			// Density average over all the water column
			rho1_aux[ind][0] += rho[zind][0]*e3[zind][0]; 
			d1_aux[ind][0]   += e3[zind][0];

			// Density average over MLD
			rho2_aux[ind][0] += (FLDARRAY)(mldmask[zind][0])*rho[zind][0]*e3[zind][0]; 
			d2_aux[ind][0]   += (FLDARRAY)(mldmask[zind][0])*e3[zind][0];			
		}

		// Second loop at the depth levels
		// Finish up computations
		for (int jj = 0; jj<mMeshPerDepth[ind].size(); ++jj) {
			int zind = mMeshPerDepth[ind][jj];

			if (cmask[zind][0] == 0) continue;

			// Barotropic Rossby radius of deformation (external)
			// Expand to all the layers
			r_ext[zind][0] = r_ext[ind][0];

			// Rossby radius of deformation (internal)
			// LR = sqrt(g'*MLD)/f_cor, where g' = g*(rho2-rho1)/rho2
			double rho2 = rho2_aux[ind][0]/d2_aux[ind][0]; // Avg over MLD
			double rho1 = rho1_aux[ind][0]/d1_aux[ind][0]; // Avg over all water column
			double gp   = this->g*fabs(rho2-rho1)/(rho1);

			r_int[zind][0] = sqrt(gp*mld[zind][0])/f_cor[ind][0]/1e3; // km
		}
	}

//	// For each depth layer, compute the mask and store MLD
//	// This algorithm requires to work on a rectilinear grid with
//	// the same number of points at each plane
//	for (int zId = 0; zId < this->zcoords.size(); zId += 1) {
//		for (int ii=0; ii<mMeshPerLayer[zId].size(); ++ii) {
//			// Indices
//			int ind0 = mMeshPerLayer[0][ii];
//			int ind  = mMeshPerLayer[zId][ii];
//			
//			if (cmask[ind][0] == 0) continue; // Mask skip
//
//			// Barotropic Rossby radius of deformation (external)
//			// LR = sqrt(g*D)/f_cor
//			r_ext[ind0][0] = sqrt(this->g*-zcoords[zId])/f_cor[ind][0]/1e3; // [km] zcoords is negative
//
//			// Density average over all the water column
//			rho1_aux[ind0][0] += rho[ind][0]*e3[ind][0]; 
//			d1_aux[ind0][0]   += e3[ind][0];
//
//			// Density average over MLD
//			rho2_aux[ind0][0] += (FLDARRAY)(mldmask[ind][0])*rho[ind][0]*e3[ind][0]; 
//			d2_aux[ind0][0]   += (FLDARRAY)(mldmask[ind][0])*e3[ind][0];
//		}
//	}
//	this->UpdateProgress(0.65);
//	
//	// Finish up computations
//	for (int zId = 0; zId < this->zcoords.size(); zId += 1) {
//		for (int ii=0; ii<mMeshPerLayer[zId].size(); ++ii) {
//			int ind0 = mMeshPerLayer[0][ii];
//			int ind  = mMeshPerLayer[zId][ii];
//
//			if (cmask[ind][0] == 0) continue;
//
//			// Barotropic Rossby radius of deformation (external)
//			// Expand to all the layers
//			r_ext[ind][0] = r_ext[ind0][0];
//
//			// Rossby radius of deformation (internal)
//			// LR = sqrt(g'*MLD)/f_cor, where g' = g*(rho2-rho1)/rho2
//			double rho2 = rho2_aux[ind0][0]/d2_aux[ind0][0]; // Avg over MLD
//			double rho1 = rho1_aux[ind0][0]/d1_aux[ind0][0]; // Avg over all water column
//			double gp   = this->g*fabs(rho2-rho1)/(rho1);
//
//			r_int[ind][0] = sqrt(gp*mld[ind][0])/f_cor[ind][0]/1e3; // km
//		}
//	}
	this->UpdateProgress(0.75);

	// Output field
	VTKARRAY *vtkr_ext, *vtkr_int;
	vtkr_ext = VTK::createVTKfromField<VTKARRAY,FLDARRAY>("Barotropic Rossby radius",r_ext);
	vtkr_int = VTK::createVTKfromField<VTKARRAY,FLDARRAY>("Internal Rossby radius",r_int);

	if (iscelld) {
		output->GetCellData()->AddArray(vtkr_ext);
		output->GetCellData()->AddArray(vtkr_int);
	} else {
		output->GetPointData()->AddArray(vtkr_ext);
		output->GetPointData()->AddArray(vtkr_int);
	}

	vtkr_ext->Delete();
	vtkr_int->Delete();

	// Return
	this->UpdateProgress(1.0);
	return 1;
}
