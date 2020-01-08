/*=========================================================================

  Program:   OGSComputeQCriterion
  Module:    vtkOGSComputeQCriterion.cxx

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkCell.h"
#include "vtkPoints.h"
#include "vtkCellData.h"
#include "vtkPointData.h"
#include "vtkFieldData.h"
#include "vtkTypeUInt8Array.h"
#include "vtkDoubleArray.h"
#include "vtkCellData.h"
#include "vtkDataArraySelection.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkRectilinearGrid.h"

#include "vtkOGSComputeQCriterion.h"

#include "vtkObjectFactory.h"

#ifdef __linux__
// Include OpenMP when working with GCC
#include <omp.h>
#define OMP_MAX_THREADS omp_get_max_threads()
#define OMP_THREAD_NUM  omp_get_thread_num()
#else
#define OMP_MAX_THREADS 1
#define OMP_THREAD_NUM  0
#endif

#ifdef PARAVIEW_USE_MPI
#include "vtkMultiProcessController.h"
vtkCxxSetObjectMacro(vtkOGSComputeQCriterion, Controller, vtkMultiProcessController);
#endif

vtkStandardNewMacro(vtkOGSComputeQCriterion);

//----------------------------------------------------------------------------

/*
	Macro to set the array precision 
*/
#define FLDARRAY double
#define FLDMASK  uint8_t
#define VTKARRAY vtkDoubleArray
#define VTKMASK  vtkTypeUInt8Array

// V3.h and field.h defined in vtkOGSDerivatives.h
#include "../_utils/fieldOperations.hpp"
#include "../_utils/vtkFields.hpp"
#include "../_utils/vtkOperations.hpp"

//----------------------------------------------------------------------------
vtkOGSComputeQCriterion::vtkOGSComputeQCriterion() {
	this->field     = NULL;
	this->coef      = 0.2;
	this->grad_type = 0;
	this->nProcs    = 0;
	this->procId    = 0;

	#ifdef PARAVIEW_USE_MPI
		this->Controller = NULL;
		this->SetController(vtkMultiProcessController::GetGlobalController());
	#endif
}

//----------------------------------------------------------------------------
vtkOGSComputeQCriterion::~vtkOGSComputeQCriterion() {
	this->Setfield(NULL);

	#ifdef PARAVIEW_USE_MPI
		this->SetController(NULL);	
	#endif
}

//----------------------------------------------------------------------------
int vtkOGSComputeQCriterion::RequestInformation(vtkInformation* vtkNotUsed(request),
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
int vtkOGSComputeQCriterion::RequestData(vtkInformation *vtkNotUsed(request),
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

	output->ShallowCopy(input);

	this->UpdateProgress(0.);

	// Recover XYZ number of nodes
	int nx = input->GetXCoordinates()->GetNumberOfTuples();
	int ny = input->GetYCoordinates()->GetNumberOfTuples();
	int nz = input->GetZCoordinates()->GetNumberOfTuples();

	// Recover a mask
	// This is useful to perform a true average and srd on the Okubo-Weiss
	// Either the basins or the coast mask will work
	field::Field<FLDMASK> mask;
	VTKMASK *vtkMask = VTKMASK::SafeDownCast( input->GetCellData()->GetArray("coast mask") );
	if (!vtkMask)
		vtkMask = VTKMASK::SafeDownCast( input->GetCellData()->GetArray("basins mask") );

	if (vtkMask)
		mask = VTK::createFieldfromVTK<VTKMASK,FLDMASK>(vtkMask);
	else
		vtkWarningMacro("Cannot find a working mask!\nResults might be affected.");

	// Recover weights
	VTKARRAY *vtke1 = VTKARRAY::SafeDownCast( input->GetCellData()->GetArray("e1") );
	VTKARRAY *vtke2 = VTKARRAY::SafeDownCast( input->GetCellData()->GetArray("e2") );
	VTKARRAY *vtke3 = VTKARRAY::SafeDownCast( input->GetCellData()->GetArray("e3") );

	if (this->grad_type > 1 && (vtke1 == NULL || vtke2 == NULL || vtke3 == NULL)) {
		vtkErrorMacro("Mesh weights (e1, e2 and e3) need to be loaded to proceed!");
		return 0;
	}

	// Convert to field arrays
	field::Field<FLDARRAY> e1, e2, e3;
	if (vtke1) e1 = VTK::createFieldfromVTK<VTKARRAY,FLDARRAY>(vtke1);
	if (vtke2) e2 = VTK::createFieldfromVTK<VTKARRAY,FLDARRAY>(vtke2);
	if (vtke3) e3 = VTK::createFieldfromVTK<VTKARRAY,FLDARRAY>(vtke3);

	// Recover velocity vector
	VTKARRAY *vtkArray;
	vtkArray = VTKARRAY::SafeDownCast(input->GetCellData()->GetArray(this->field));
	if (vtkArray == NULL) {
		vtkErrorMacro("Input array "<<this->field<<"is not a cell array."<<
			"This function can only deal with cell arrays!");
		return 0;
	}

	// Here we have a valid vtkArray. Convert it to a field.
	field::Field<FLDARRAY> array;
	array = VTK::createFieldfromVTK<VTKARRAY,FLDARRAY>(vtkArray);

	// Cell centers, update only when request info
	if (this->isReqInfo) {
		this->isReqInfo = false;
		
		// Recover Metadata array (depth factor)
		vtkStringArray *vtkmetadata = vtkStringArray::SafeDownCast(
			input->GetFieldData()->GetAbstractArray("Metadata"));
		double dfact = (vtkmetadata != NULL) ? std::stod( vtkmetadata->GetValue(2) ) : 1000.;
		
		if (vtkmetadata == NULL) 
			vtkWarningMacro("Field array Metadata not found! Depth factor set to 1000. automatically.");

		this->xyz = VTK::getVTKCellCenters(input,dfact);
	}

	// Loop on the 3D mesh and compute the Q criterion also new scalar arrays 
	// to store the values for the weights, the mean and thestandard deviation.
	// This code uses the recursive updating formulae from:
	// Chan, Tony F., Gene H. Golub, and Randall J. LeVeque. "Updating formulae and a pairwise algorithm for computing sample variances." 
	// COMPSTAT 1982 5th Symposium held at Toulouse 1982. Physica, Heidelberg, 1982.
	// To compute the statistics. For a bit of memory and maybe a fairly slower first loop
	// we are able to avoid doing a second mesh loop.
	field::Field<FLDARRAY> Q(array.get_n(),1,0.), Q_stats(3,OMP_MAX_THREADS,0.);

	this->UpdateProgress(0.1);

	#pragma omp parallel for collapse(3)
	for (int kk = 0; kk < nz-1; kk++) {
		for (int jj = 0; jj < ny-1; jj++) {
			for (int ii = 0; ii < nx-1; ii++) {
				// Compute current point
				int ind = CLLIND(ii,jj,kk,nx,ny);

				// Use the mask to determine whether we are in the sea or not
				// to improve the averaging
				if (!mask.isempty()) {
					bool skip_ind = true;
					for (int ii = 0; ii < mask.get_m(); ++ii) {
						if (mask[ind][ii] > 0) {skip_ind = false; break;}
					}
					if (skip_ind) continue;
				}

				// Compute the gradient
				FLDARRAY deri[9] = {0.,0.,0.,0.,0.,0.,0.,0.,0.};

				switch (this->grad_type) {
					case 0: // Second order, face centered gradient
							// This gradient is unsafe as it relies on the mesh projection
						field::gradXYZ2_ijk(ii,jj,kk,nx-1,ny-1,nz-1,this->xyz,array,deri);
						break;
					case 1: // Fourth order, face centered gradient
							// This gradient is unsafe as it relies on the mesh projection
						field::gradXYZ4_ijk(ii,jj,kk,nx-1,ny-1,nz-1,this->xyz,array,deri);
						break;
					case 2:	// OGSTM-BFM approach according to the NEMO handbook
							// This gradient is safe as it relies on the code implementation
						field::gradOGS1_ijk(ii,jj,kk,nx-1,ny-1,nz-1,array,e1,e2,e3,deri);
						break;
					case 3:	// 2nd order OGSTM-BFM approach
							// This gradient is experimental
						field::gradOGS2_ijk(ii,jj,kk,nx-1,ny-1,nz-1,array,e1,e2,e3,deri);
						break;
					case 4:	// 4th order OGSTM-BFM approach
							// This gradient is experimental
						field::gradOGS4_ijk(ii,jj,kk,nx-1,ny-1,nz-1,array,e1,e2,e3,deri);
						break;
				}

				// Compute Q-criterion
				Q[ind][0] = -0.5*(deri[0]*deri[0] + deri[4]*deri[4] + deri[8]*deri[8])
						    -deri[1]*deri[3]-deri[2]*deri[6]-deri[5]*deri[7];
				Q[ind][0]*= -2.; // Scaling fix so that Q == OW

				// Compute the statistics
				int        tId = OMP_THREAD_NUM;  // Which thread are we
				FLDARRAY     w = e1[ind][0]*e2[ind][0]; // Weight
				FLDARRAY m_old = Q_stats[1][tId];      // Mean Old

				Q_stats[0][tId] += w;                                                   // sum(weight)
				Q_stats[1][tId] += w/Q_stats[0][tId]*(Q[ind][0]-Q_stats[1][tId]);       // accumulate mean
				Q_stats[2][tId] += w*(Q[ind][0] - m_old)*(Q[ind][0] - Q_stats[1][tId]); // accumulate std deviation
			}
		}
	}

	// Q_stats is a shared array containing the statistics per each thread. We must now
	// reduce to obtain the mean and standard deviation
	FLDARRAY sum_weights = 0., Q_mean = 0., Q_std = 0.;
	for (int ii = 0; ii < OMP_MAX_THREADS; ++ii) {
		FLDARRAY sum_weights_old = sum_weights;
		// Weights
		sum_weights += Q_stats[0][ii];
		// Mean
		Q_mean = (sum_weights_old*Q_mean + Q_stats[0][ii]*Q_stats[1][ii])/sum_weights;
		// Standard deviation
		FLDARRAY delta = Q_stats[1][ii] - Q_mean;
		Q_std += Q_stats[2][ii] + delta*delta*sum_weights_old*Q_stats[0][ii]/sum_weights;
	}
	Q_std = this->coef*sqrt(Q_std/sum_weights);

	this->UpdateProgress(0.4);

	// We now have the standard deviation, let us work out the mask
	field::Field<FLDMASK> Qm(array.get_n(),1);

	#pragma omp parallel for collapse(3)
	for (int kk = 0; kk < nz-1; kk++) {
		for(int jj = 0; jj < ny-1; jj++){
			for (int ii = 0; ii < nx-1; ii++) {
				// Compute current point
				int ind = CLLIND(ii,jj,kk,nx,ny);

				// Initalize mask for Background field
				Qm[ind][0] = 1;

				// Use the mask to determine whether we are in the sea or not
				// to improve the averaging
				if (!mask.isempty()) {
					bool skip_ind = true;
					for (int ii = 0; ii < mask.get_m(); ++ii) {
						if (mask[ind][ii] > 0) {skip_ind = false; break;}
					}
					if (skip_ind) continue;
				}

				// Set the values of the mask
				if (Q[ind][0] < -Q_std) Qm[ind][0] = 0; // Vorticity-dominated flow
				if (Q[ind][0] >  Q_std) Qm[ind][0] = 2; // Strain-dominated flow

				// Undo the scaling to recover the definition of Q
				Q[ind][0]*= -.5;
			}
		}
	}
	this->UpdateProgress(0.9);
		
	// Generate VTK arrays and add them to the output
	vtkArray = VTK::createVTKfromField<VTKARRAY,FLDARRAY>("Qcriterion",Q);
	output->GetCellData()->AddArray(vtkArray);  vtkArray->Delete();
	vtkMask = VTK::createVTKfromField<VTKMASK,FLDMASK>("Qcriterion_mask",Qm);
	output->GetCellData()->AddArray(vtkMask);  vtkMask->Delete();

	// Make "Qcriterion" as the active scalar
	output->GetCellData()->SetActiveScalars("Qcriterion");

	// Copy the input grid
	this->UpdateProgress(1.);
	return 1;
}
