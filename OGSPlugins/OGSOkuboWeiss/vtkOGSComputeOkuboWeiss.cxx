/*=========================================================================

  Program:   OGSComputeOkuboWeiss
  Module:    vtkOGSComputeOkuboWeiss.cxx

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkOGSComputeOkuboWeiss.h"

#include "vtkCell.h"
#include "vtkPoints.h"
#include "vtkCellData.h"
#include "vtkPointData.h"
#include "vtkFieldData.h"
#include "vtkTypeUInt8Array.h"
#include "vtkFloatArray.h"
#include "vtkDoubleArray.h"
#include "vtkCellData.h"
#include "vtkDataArraySelection.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkRectilinearGrid.h"

#include "vtkObjectFactory.h"

#include <cstdint>
#include <omp.h>
int omp_get_max_threads();
int omp_get_thread_num();


#ifdef PARAVIEW_USE_MPI
#include "vtkMultiProcessController.h"
vtkCxxSetObjectMacro(vtkOGSComputeOkuboWeiss, Controller, vtkMultiProcessController);
#endif

vtkStandardNewMacro(vtkOGSComputeOkuboWeiss);

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
vtkOGSComputeOkuboWeiss::vtkOGSComputeOkuboWeiss() {
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
vtkOGSComputeOkuboWeiss::~vtkOGSComputeOkuboWeiss() {
	this->Setfield(NULL);

	#ifdef PARAVIEW_USE_MPI
		this->SetController(NULL);	
	#endif
}

//----------------------------------------------------------------------------
int vtkOGSComputeOkuboWeiss::RequestInformation(vtkInformation* vtkNotUsed(request),
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
int vtkOGSComputeOkuboWeiss::RequestData(vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector, vtkInformationVector *outputVector) {
	// Get the info objects
	vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
	vtkInformation *outInfo = outputVector->GetInformationObject(0);

	// Get the input and output
	vtkRectilinearGrid *input = vtkRectilinearGrid::SafeDownCast(
		inInfo->Get(vtkDataObject::DATA_OBJECT()));
	vtkRectilinearGrid *output = vtkRectilinearGrid::SafeDownCast(
		outInfo->Get(vtkDataObject::DATA_OBJECT()));

	// Stop all threads except from the master to execute
	#ifdef PARAVIEW_USE_MPI
	if (this->procId > 0) return 1;
	#endif

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

	this->UpdateProgress(0.1);

	// Loop on the 2D surface mesh (k == 0) and compute the Okubo-Weiss parameter
	// also new scalar arrays to store the Okubo-Weiss values for the weights, mean
	// and standard deviation.
	// This code uses the recursive updating formulae from:
	// Chan, Tony F., Gene H. Golub, and Randall J. LeVeque. "Updating formulae and a pairwise algorithm for computing sample variances." 
	// COMPSTAT 1982 5th Symposium held at Toulouse 1982. Physica, Heidelberg, 1982.
	// To compute the statistics. For a bit of memory and maybe a fairly slower first loop
	// we are able to avoid doing a second mesh loop.
	field::Field<FLDARRAY> OW(array.get_n(),1,0.), OW_stats(3,omp_get_max_threads(),0.);

	#pragma omp parallel for collapse(2)
	for (int jj = 0; jj < ny; ++jj) { 
		for (int ii = 0; ii < nx; ++ii) {
			// Compute current point
			int ind = CLLIND(ii,jj,0,nx,ny);

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
					field::gradXYZ2_ijk(ii,jj,0,nx-1,ny-1,nz-1,this->xyz,array,deri);
					break;
				case 1: // Fourth order, face centered gradient
						// This gradient is unsafe as it relies on the mesh projection
					field::gradXYZ4_ijk(ii,jj,0,nx-1,ny-1,nz-1,this->xyz,array,deri);
					break;
				case 2:	// OGSTM-BFM approach according to the NEMO handbook
						// This gradient is safe as it relies on the code implementation
					field::gradOGS1_ijk(ii,jj,0,nx-1,ny-1,nz-1,array,e1,e2,e3,deri);
					break;
				case 3:	// 2nd order OGSTM-BFM approach
						// This gradient is experimental
					field::gradOGS2_ijk(ii,jj,0,nx-1,ny-1,nz-1,array,e1,e2,e3,deri);
					break;
				case 4:	// 4th order OGSTM-BFM approach
						// This gradient is experimental
					field::gradOGS4_ijk(ii,jj,0,nx-1,ny-1,nz-1,array,e1,e2,e3,deri);
					break;
			}

			// Rates of strain and stress
			FLDARRAY Sn = (deri[0] - deri[4]); // dudx - dvdy
			FLDARRAY Ss = (deri[3] + deri[1]); // dvdx + dudy
			FLDARRAY W  = (deri[3] - deri[1]); // dvdx - dudy

			// Compute Okubo-Weiss
			OW[ind][0] = Sn*Sn + Ss*Ss - W*W;
			
			// Compute the statistics
			int        tId = omp_get_thread_num();  // Which thread are we
			FLDARRAY     w = e1[ind][0]*e2[ind][0]; // Weight
			FLDARRAY m_old = OW_stats[1][tId];      // Mean Old

			OW_stats[0][tId] += w;                                                      // sum(weight)
			OW_stats[1][tId] += w/OW_stats[0][tId]*(OW[ind][0]-OW_stats[1][tId]);       // accumulate mean
			OW_stats[2][tId] += w*(OW[ind][0] - m_old)*(OW[ind][0] - OW_stats[1][tId]); // accumulate std deviation
		}
	}
	this->UpdateProgress(0.2);

	// OW_stats is a shared array containing the statistics per each thread. We must now
	// reduce to obtain the mean and standard deviation
	FLDARRAY sum_weights = 0., OW_mean = 0., OW_std = 0.;
	for (int ii = 0; ii < omp_get_max_threads(); ++ii) {
		FLDARRAY sum_weights_old = sum_weights;
		// Weights
		sum_weights += OW_stats[0][ii];
		// Mean
		OW_mean = (sum_weights_old*OW_mean + OW_stats[0][ii]*OW_stats[1][ii])/sum_weights;
		// Standard deviation
		FLDARRAY delta = OW_stats[1][ii] - OW_mean;
		OW_std += OW_stats[2][ii] + delta*delta*sum_weights_old*OW_stats[0][ii]/sum_weights;
	}
	OW_std = this->coef*sqrt(OW_std/sum_weights);

	this->UpdateProgress(0.4);

	// We now have the standard deviation, let us do a full mesh loop to set the 
	// values of the mask and the OW array
	field::Field<FLDMASK> OWm(array.get_n(),1);

	#pragma omp parallel for collapse(3) 
	for (int kk = 0; kk < nz-1; kk++) {
		for(int jj = 0; jj < ny-1; jj++) {
			for (int ii = 0; ii < nx-1; ii++) {
				// Compute current point
				int ind0 = CLLIND(ii,jj,0,nx,ny);
				int ind  = CLLIND(ii,jj,kk,nx,ny);

				// Initalize mask for Background field
				OWm[ind][0] = 1;

				// Use the mask to determine whether we are in the sea or not.
				// This will make the arrays follow the batimetry
				if (!mask.isempty()) {
					bool skip_ind = true;
					for (int ii = 0; ii < mask.get_m(); ++ii) {
						if (mask[ind][ii] > 0) {skip_ind = false; break;}
					}
					if (skip_ind) continue;
				}

				// Set the Okubo-Weiss array
				OW[ind][0]  = OW[ind0][0];

				// Set the values of the mask
				if (OW[ind0][0] < -OW_std) OWm[ind][0] = 0; // Vorticity-dominated flow
				if (OW[ind0][0] >  OW_std) OWm[ind][0] = 2; // Strain-dominated flow
			}
		}
	}

	this->UpdateProgress(0.8);

	// Generate VTK arrays and add them to the output
	vtkArray = VTK::createVTKfromField<VTKARRAY,FLDARRAY>("Okubo-Weiss",OW);
	output->GetCellData()->AddArray(vtkArray);  vtkArray->Delete();
	vtkMask = VTK::createVTKfromField<VTKMASK,FLDMASK>("Okubo-Weiss mask",OWm);
	output->GetCellData()->AddArray(vtkMask);  vtkMask->Delete();

	// Copy the input grid
	this->UpdateProgress(1.);
	return 1;
}
