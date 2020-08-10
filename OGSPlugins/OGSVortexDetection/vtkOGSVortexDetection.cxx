/*=========================================================================

  Program:   OGSVortexDetection
  Module:    vtkOGSVortexDetection.cxx

  Copyright (c) 2020 Arnau Miro, OGS
  All rights reserved.

	 This software is distributed WITHOUT ANY WARRANTY; without even
	 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
	 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkOGSVortexDetection.h"

#include "vtkIdList.h"
#include "vtkCell.h"
#include "vtkPoints.h"
#include "vtkCellData.h"
#include "vtkPointData.h"
#include "vtkFieldData.h"
#include "vtkTypeUInt8Array.h"
#include "vtkFloatArray.h"
#include "vtkDoubleArray.h"
#include "vtkTable.h"
#include "vtkDataArraySelection.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#include "vtkObjectFactory.h"

#include <vector>

#ifdef PARAVIEW_USE_MPI
#include "vtkMultiProcessController.h"
vtkCxxSetObjectMacro(vtkOGSVortexDetection, Controller, vtkMultiProcessController);
#endif

vtkStandardNewMacro(vtkOGSVortexDetection);

//----------------------------------------------------------------------------
#include "macros.h"
#include "field.h"
#include "projection.h"
#include "vtkFields.h"
#include "vtkOperations.h"

//----------------------------------------------------------------------------
typedef enum _MASK_TYPE {
	NONE,
	OKUBOWEISS,
	QCRITERION,
	RORTEX
}MASK_TYPE;

MASK_TYPE int2mask(const int imask) {
	if (imask == 0) return OKUBOWEISS;
	if (imask == 1) return QCRITERION;
	if (imask == 2) return RORTEX;
	return NONE;
}

MASK_TYPE str2mask(const char *mask) {
	std::string strmask(mask);
	std::transform(strmask.begin(), strmask.end(), strmask.begin(), ::tolower);
	if (strmask == "qcriterion")  return QCRITERION;
	if (strmask == "q-criterion") return QCRITERION;
	if (strmask == "okuboweiss")  return OKUBOWEISS;
	if (strmask == "okubo-weiss") return OKUBOWEISS;
	if (strmask == "rortex")      return RORTEX;
	return NONE;
}

//----------------------------------------------------------------------------
vtkOGSVortexDetection::vtkOGSVortexDetection() {
	this->SetNumberOfInputPorts(1);
	this->SetNumberOfOutputPorts(2);

	this->vecFName     = NULL;
	this->scaFName     = NULL;
//	this->normName     = NULL;
	this->arr2Mask     = NULL;
	this->maskName     = NULL;

	this->projName     = std::string("Mercator");
	this->dfact        = 1000.;
	this->computeMask  = 1;
	this->coef         = 0.2;
	this->maxreps      = 5;
	this->maxiter      = 100;
	this->minres       = 5;
	this->changemask   = 1;

	this->computeBCent = 0;
	this->computeCent  = 1;
	this->computeSize  = 1;
	this->computeRAxis = 1;
	this->computeAStr  = 1;
	this->computeRStr  = 1;
	this->computeCirc  = 0;

	this->nProcs       = 0;
	this->procId       = 0;

	this->Projections  = vtkStringArray::New();
	
	#ifdef PARAVIEW_USE_MPI
		this->Controller = NULL;
		this->SetController(vtkMultiProcessController::GetGlobalController());
	#endif
}

//----------------------------------------------------------------------------
vtkOGSVortexDetection::~vtkOGSVortexDetection() {
	this->SetvecFName(NULL);
	this->SetscaFName(NULL);
//	this->SetnormName(NULL);
	this->Setarr2Mask(NULL);
	this->SetmaskName(NULL);

	#ifdef PARAVIEW_USE_MPI
		this->SetController(NULL);	
	#endif	
}

//----------------------------------------------------------------------------
int vtkOGSVortexDetection::RequestInformation(vtkInformation* vtkNotUsed(request),
	vtkInformationVector** vtkNotUsed(inputVector), vtkInformationVector* vtkNotUsed(outputVector)) {

	#ifdef PARAVIEW_USE_MPI
	if (this->Controller->GetNumberOfProcesses() > 1) {
		this->nProcs = this->Controller->GetNumberOfProcesses();
		this->procId = this->Controller->GetLocalProcessId();
	}

	// Stop all threads except from the master to execute
	if (this->procId > 0) return 1;
	#endif

	return 1;
}

//----------------------------------------------------------------------------
int vtkOGSVortexDetection::RequestData(vtkInformation* vtkNotUsed(request),
	vtkInformationVector** inputVector, vtkInformationVector* outputVector ) {
	// get the info objects
	vtkInformation *inInfo   = inputVector[0]->GetInformationObject(0);
	vtkInformation *outInfoA = outputVector->GetInformationObject(0);
	vtkInformation *outInfoB = outputVector->GetInformationObject(1);

	// Stop all threads except from the master to execute
	#ifdef PARAVIEW_USE_MPI
	if (this->procId > 0) return 1;
	#endif

	// get the input and output
	vtkDataSet *input = vtkDataSet::SafeDownCast(
		inInfo->Get(vtkDataObject::DATA_OBJECT()));
	vtkDataSet *outputA = vtkDataSet::SafeDownCast(
		outInfoA->Get(vtkDataObject::DATA_OBJECT()));
	vtkTable *outputB = vtkTable::SafeDownCast(
		outInfoB->Get(vtkDataObject::DATA_OBJECT()));
	
	outputA->ShallowCopy(input);

	// Decide whether we have cell or point data
	int n_cell_vars  = outputA->GetCellData()->GetNumberOfArrays();
	int n_point_vars = outputA->GetPointData()->GetNumberOfArrays();

	bool iscelld = (n_cell_vars > n_point_vars) ? true : false;

	// Compute the mask field if requested
	if (this->computeMask) {
		if (this->computeMaskArray(iscelld,outputA) == 0) 
			return 0;
	}
	this->UpdateProgress(0.4);

	// Compute the vortex detection field
	vortex::VortexList vortexList;
	int nvortices = this->computeDetectionMask(iscelld,outputA,vortexList);
	if (nvortices == 0) return 0;
	this->UpdateProgress(0.7);

	// Compute vortex properties
	if (this->computeVortexProperties(iscelld,vortexList,outputA,outputB) == 0) 
		return 0;
	outputB->GetFieldData()->RemoveArray("Metadata");

	// Update progress and leave
	this->UpdateProgress(1.);
	return 1;
}

//----------------------------------------------------------------------------
vtkDataSet* vtkOGSVortexDetection::GetOutputA() {
	return vtkDataSet::SafeDownCast(this->GetOutputDataObject(0));
}

//----------------------------------------------------------------------------
vtkTable* vtkOGSVortexDetection::GetOutputB() {
	return vtkTable::SafeDownCast(this->GetOutputDataObject(1));
}

//----------------------------------------------------------------------------
void vtkOGSVortexDetection::SetOutput(vtkDataObject* d) {
	this->GetExecutive()->SetOutputData(0,d);
}

//----------------------------------------------------------------------------
int vtkOGSVortexDetection::FillOutputPortInformation(int port, vtkInformation* info) {

	if (port == 0)
		info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataSet");
	else if (port == 1)
		info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkTable");

	return 1;
}

//----------------------------------------------------------------------------
int vtkOGSVortexDetection::RequestDataObject(vtkInformation* vtkNotUsed(request),
	vtkInformationVector** inputVector, vtkInformationVector* outputVector ) {

	vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
	if (!inInfo)
		return 0;

	vtkDataSet* input = vtkDataSet::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
	if (input) {
		// Port 0
		vtkInformation* outInfoA = outputVector->GetInformationObject(0);
		vtkDataSet*      outputA = vtkDataSet::SafeDownCast(outInfoA->Get(vtkDataObject::DATA_OBJECT()));

		if (!outputA) {
			vtkDataSet* newOutputA = input->NewInstance();
			outInfoA->Set(vtkDataObject::DATA_OBJECT(), newOutputA);
			newOutputA->Delete();
		}

		// Port 1
		vtkInformation* outInfoB = outputVector->GetInformationObject(1);
		vtkTable*        outputB = vtkTable::SafeDownCast(outInfoB->Get(vtkDataObject::DATA_OBJECT()));

		if (!outputB) {
			vtkTable* newOutputB = vtkTable::New();
			outInfoB->Set(vtkDataObject::DATA_OBJECT(), newOutputB);
			newOutputB->Delete();
		}

		return 1;
	}

	return 0;
}

//----------------------------------------------------------------------------
int vtkOGSVortexDetection::RequestUpdateExtent(vtkInformation* vtkNotUsed(request),
	vtkInformationVector** inputVector, vtkInformationVector* vtkNotUsed(outputVector)) {

	int numInputPorts = this->GetNumberOfInputPorts();
	for (int i=0; i<numInputPorts; i++) {
		int numInputConnections = this->GetNumberOfInputConnections(i);
		for (int j=0; j<numInputConnections; j++) {
			vtkInformation* inputInfo = inputVector[i]->GetInformationObject(j);
			inputInfo->Set(vtkStreamingDemandDrivenPipeline::EXACT_EXTENT(), 1);
		}
	}

	return 1;
}

//----------------------------------------------------------------------------
int vtkOGSVortexDetection::computeMaskArray(bool iscelld, vtkDataSet *d) {
	
	// Load array as cell or point data
	VTKARRAY *vtkArray;
	if (iscelld)
		vtkArray = VTKARRAY::SafeDownCast( d->GetCellData()->GetArray(this->arr2Mask) );
	else
		vtkArray = VTKARRAY::SafeDownCast( d->GetPointData()->GetArray(this->arr2Mask) );

	if (vtkArray == NULL) {
		vtkErrorMacro("Cannot load mask field "<<this->arr2Mask<<"!");
		return 0;
	}

	field::Field<FLDARRAY> array = VTK::createFieldfromVTK<VTKARRAY,FLDARRAY>(vtkArray);

	// Recover weights
	VTKARRAY *vtke1, *vtke2, *vtke3;
	if (iscelld) {
		vtke1 = VTKARRAY::SafeDownCast( d->GetCellData()->GetArray("e1") );
		vtke2 = VTKARRAY::SafeDownCast( d->GetCellData()->GetArray("e2") );
		vtke3 = VTKARRAY::SafeDownCast( d->GetCellData()->GetArray("e3") );
	} else {
		vtke1 = VTKARRAY::SafeDownCast( d->GetPointData()->GetArray("e1") );
		vtke2 = VTKARRAY::SafeDownCast( d->GetPointData()->GetArray("e2") );
		vtke3 = VTKARRAY::SafeDownCast( d->GetPointData()->GetArray("e3") );		
	}

	if (vtke1 == NULL || vtke2 == NULL || vtke3 == NULL) {
		vtkErrorMacro("Mesh weights (e1, e2 and e3) need to be loaded to proceed!");
		return 0;
	}

	// Convert to field arrays
	field::Field<FLDARRAY> e1, e2, e3;
	if (vtke1) e1 = VTK::createFieldfromVTK<VTKARRAY,FLDARRAY>(vtke1);
	if (vtke2) e2 = VTK::createFieldfromVTK<VTKARRAY,FLDARRAY>(vtke2);
	if (vtke3) e3 = VTK::createFieldfromVTK<VTKARRAY,FLDARRAY>(vtke3);

	// Loop the mesh and compute the standard deviation
	field::Field<FLDARRAY> stats(3,OMP_MAX_THREADS,0.);
	#pragma omp parallel
	{
	for(int ii=OMP_THREAD_NUM; ii<array.get_n(); ii+=OMP_NUM_THREADS) {

		FLDARRAY val = (array.get_m() == 1) ? array[ii][0] : 
			sqrt(array[ii][0]*array[ii][0] + array[ii][1]*array[ii][1] + array[ii][2]*array[ii][2]);

		// Compute the statistics
		int        tId = OMP_THREAD_NUM;      // Which thread are we
		FLDARRAY     w = e1[ii][0]*e2[ii][0]; // Weight
		FLDARRAY m_old = stats[1][tId];       // Mean Old

		stats[0][tId] += w;                                     // sum(weight)
		stats[1][tId] += w/stats[0][tId]*(val - stats[1][tId]); // accumulate mean
		stats[2][tId] += w*(val - m_old)*(val - stats[1][tId]); // accumulate std deviation
	}
	} // end omp
	this->UpdateProgress(0.1);

	// stats is a shared array containing the statistics per each thread. We must now
	// reduce to obtain the mean and standard deviation
	FLDARRAY sum_weights = 0., mean = 0., std = 0.;
	for (int ii = 0; ii < OMP_MAX_THREADS; ++ii) {
		FLDARRAY sum_weights_old = sum_weights;
		// Weights
		sum_weights += stats[0][ii];
		// Mean
		mean = (sum_weights_old*mean + stats[0][ii]*stats[1][ii])/sum_weights;
		// Standard deviation
		FLDARRAY delta = stats[1][ii] - mean;
		std += stats[2][ii] + delta*delta*sum_weights_old*stats[0][ii]/sum_weights;
	}
	std = this->coef*sqrt(std/sum_weights);
	this->UpdateProgress(0.2);

	// Second loop on the mesh to set the mask values
	field::Field<FLDMASK> mask(array.get_n(),1,(FLDMASK)(0));
	#pragma omp parallel
	{
	for(int ii=OMP_THREAD_NUM; ii<array.get_n(); ii+=OMP_NUM_THREADS) {

		FLDARRAY val = (array.get_m() == 1) ? array[ii][0] : 
			sqrt(array[ii][0]*array[ii][0] + array[ii][1]*array[ii][1] + array[ii][2]*array[ii][2]);

		// Set the values of the mask
		switch (str2mask(this->arr2Mask)) {
			case OKUBOWEISS:
				mask[ii][0] = (val < -std) ? (FLDMASK)(1) : 0; // Vorticity-dominated flow
				break;
			case QCRITERION:
				mask[ii][0] = (val > 0.5*std) ? (FLDMASK)(1) : 0; // Vorticity-dominated flow
				break;
			case RORTEX:
				mask[ii][0] = (val > std)     ? (FLDMASK)(1) : 0; // Vorticity-dominated flow
				break;
			default:
				mask[ii][0] = 0;
				break;
		}
	}
	} // end omp
	this->UpdateProgress(0.3);

	// Add mask to output
	this->SetmaskName("mask");
	VTKMASK *vtkMask; vtkMask = VTK::createVTKfromField<VTKMASK,FLDMASK>(this->maskName,mask);

	if (iscelld)
		d->GetCellData()->AddArray(vtkMask);
	else
		d->GetPointData()->AddArray(vtkMask);

	vtkMask->Delete();

	return 1;
}

//----------------------------------------------------------------------------
int vtkOGSVortexDetection::computeDetectionMask(bool iscelld, vtkDataSet *d, vortex::VortexList &vortexList) {
	int nvortices = 0;
	// Load mask as cell or point data
	VTKMASK *vtkMask;
	if (iscelld)
		vtkMask = VTKMASK::SafeDownCast( d->GetCellData()->GetArray(this->maskName) );
	else
		vtkMask = VTKMASK::SafeDownCast( d->GetPointData()->GetArray(this->maskName) );

	if (!vtkMask) {
		vtkErrorMacro("Error loading <"<<this->maskName<<">! Aborting...");
		return 0;
	}

	field::Field<FLDMASK> mask = VTK::createFieldfromVTK<VTKMASK,FLDMASK>(vtkMask);
	field::Field<FLDARRAY> lmask(mask.get_n(),1,0.);

	// Loop the mesh 
	int clabel = 0;
	for (int ii = 0; ii < mask.get_n(); ++ii) {
		int plabel = 0;
		// Discard the point if we are not in a marked zone
		if (!mask[ii][0]) continue;
		// If the point is valid check its neighbours
		std::vector<int> neighbours = (iscelld) ? VTK::getConnectedCells(d,ii) : VTK::getConnectedPoints(d,ii);
		// Loop the neighbours
		for (int jj : neighbours) {
			// Discard invalid point
			if (!mask[jj][0]) continue;
			// Set the paint label
			plabel = (lmask[jj][0] > 0. && lmask[jj][0] < lmask[ii][0]) ? (int)(lmask[jj][0]) : plabel;
		}
		// If we haven't found any suitable label increase the label
		if (!plabel) { clabel++; plabel = clabel; }
		// Paint the current point and its neighbours
		lmask[ii][0] = (FLDARRAY)(plabel);
		for (int jj : neighbours) {
			if (!mask[jj][0]) continue;
			lmask[jj][0] = (FLDARRAY)(plabel);
		}
	}
	this->UpdateProgress(0.5);

	nvortices = clabel;

	// Second loop
	int iter = 0, reps = 0;
	for (iter=0; iter<this->maxiter && reps<this->maxreps; ++iter) {
		int nvortices_old = nvortices;
		nvortices = 0;
		for (int ii = 0; ii < mask.get_n(); ++ii) {
			// Reverse every iteration
			int ind = (iter%2 == 0) ? (mask.get_n()-1) - ii : ii;
			// Discard the point if we are not in a marked zone
			if (!mask[ind][0]) continue;
			// If the point is valid check its neighbours
			std::vector<int> neighbours = (iscelld) ? VTK::getConnectedCells(d,ind) : VTK::getConnectedPoints(d,ind);
			// Loop the neighbours
			int plabel = (int)(lmask[ind][0]);
			for (int jj : neighbours) {
				// Discard invalid point
				if (!mask[jj][0]) continue;
				// Set the paint label
				plabel = (lmask[jj][0] > 0. && lmask[jj][0] < lmask[ind][0]) ? (int)(lmask[jj][0]) : plabel;
			}
			// Paint the current point and its neighbours
			lmask[ind][0] = (FLDARRAY)(plabel);
			nvortices    = (plabel > nvortices) ? plabel : nvortices;
			for (int jj : neighbours) {
				if (!mask[jj][0]) continue;
				lmask[jj][0] = (FLDARRAY)(plabel);
			}
		}
		// Exit criteria
		if (nvortices == nvortices_old)
			reps++;
		else
			reps = 0;
	}

	// Stop if we didn't find any vortex
	if (nvortices == 0)
		return 0;
	else
		nvortices++; // Since C starts at 0

	this->UpdateProgress(0.6);

	// Loop the mesh yet again, this time we will generate a vector for each vortex and
	// we will save the elements that belong to said vortex
	std::vector<std::vector<int>> vortex_elems(nvortices);

	for (int ii = 0; ii < mask.get_n(); ++ii) {
		// Discard the point if we are not in a marked zone
		if (!mask[ii][0]) continue;
		// If the point is not discarded it means we are in a vortex,
		// recover the vortex id
		int ivort = (int)(lmask[ii][0]);
		// Store the element id
		vortex_elems[ivort].push_back(ii);
	}

	// Loop on the number of vortices and delete those that do not fit
	// the minimum criteria to be a vortex
	for (int ivort = 0; ivort < nvortices; ++ivort) {
		int nvort = (int)(vortex_elems[ivort].size());
		// Decide whether we have a vortex with enough resolution
		// or not
		if (nvort <= this->minres) {
			// Vortex does not fit the minimum resolution criteria
			// erase it from lmask (mask keeps all the values)
			for (int iel : vortex_elems[ivort])
				lmask[iel][0] = 0;
		} else {
			// Set up the vortex structure
			vortexList.push_back( 
				vortex::Vortex((int)(vortex_elems[ivort].size()),vortex_elems[ivort].data()) 
				);
		}
	}

	// Number of vortices
	nvortices = (int)(vortexList.size());

	// Re-enumerate the existing vortices
	for (int ivort = 0; ivort < nvortices; ++ivort) {
		// Loop all the elements that form the vortex
		for (int ii = 0; ii < vortexList[ivort].get_n(); ++ii) {
			int iel = vortexList[ivort][ii];
			lmask[iel][0] = ivort+1; // Since we start counting at 0
		}
	}

	// Add mask to output
	VTKARRAY *vtkArray;
	vtkArray = VTK::createVTKfromField<VTKARRAY,FLDARRAY>("vortex mask",lmask);

	if (iscelld)
		d->GetCellData()->AddArray(vtkArray);
	else
		d->GetPointData()->AddArray(vtkArray);

	vtkArray->Delete();
	return nvortices;
}

//----------------------------------------------------------------------------
int vtkOGSVortexDetection::computeVortexProperties(bool iscelld, vortex::VortexList &vortexList, 
	vtkDataSet *d, vtkTable *t) {

	// Recover Metadata array (depth factor)
	vtkStringArray *vtkmetadata = vtkStringArray::SafeDownCast(
		d->GetFieldData()->GetAbstractArray("Metadata"));
	
	if (vtkmetadata == NULL)
		vtkWarningMacro("Field array Metadata not found! Using user input values.");
	
	this->dfact    = (vtkmetadata != NULL) ? std::stod( vtkmetadata->GetValue(2) ) : this->dfact;
	this->projName = (vtkmetadata != NULL) ? vtkmetadata->GetValue(7) : std::string("Mercator");
	std::transform(this->projName.begin(), this->projName.end(), this->projName.begin(), ::tolower);

	v3::V3v xyz = (iscelld) ? VTK::getVTKCellCenters(d,this->dfact) : VTK::getVTKCellPoints(d,this->dfact);

	// Load arrays
	VTKMASK *vtkMask;
	VTKARRAY *vtkLmsk, *vtkScaf, *vtkVecf, *vtkColumn;
	if (iscelld) {
		vtkScaf = VTKARRAY::SafeDownCast( d->GetCellData()->GetArray(this->scaFName) );
		vtkVecf = VTKARRAY::SafeDownCast( d->GetCellData()->GetArray(this->vecFName) );
		vtkLmsk = VTKARRAY::SafeDownCast( d->GetCellData()->GetArray("vortex mask") );
		vtkMask = VTKMASK::SafeDownCast(  d->GetCellData()->GetArray(this->maskName) );
	} else {
		vtkScaf = VTKARRAY::SafeDownCast( d->GetPointData()->GetArray(this->scaFName) );
		vtkVecf = VTKARRAY::SafeDownCast( d->GetPointData()->GetArray(this->vecFName) );
		vtkLmsk = VTKARRAY::SafeDownCast( d->GetPointData()->GetArray("vortex mask") );
		vtkMask = VTKMASK::SafeDownCast(  d->GetPointData()->GetArray(this->maskName) );
	}

	if (vtkScaf == NULL) {
		vtkErrorMacro("Cannot load field <"<<this->scaFName<<">! Aborting...");
		return 0;
	}
	if (vtkVecf == NULL) {
		vtkErrorMacro("Cannot load field <"<<this->vecFName<<">! Aborting...");
		return 0;
	}
	if (vtkMask == NULL) {
		vtkErrorMacro("Cannot load field <"<<this->maskName<<">! Aborting...");
		return 0;
	}

	field::Field<FLDARRAY> scaf  = VTK::createFieldfromVTK<VTKARRAY,FLDARRAY>(vtkScaf);
	field::Field<FLDARRAY> vecf  = VTK::createFieldfromVTK<VTKARRAY,FLDARRAY>(vtkVecf);
	field::Field<FLDARRAY> lmask = VTK::createFieldfromVTK<VTKARRAY,FLDARRAY>(vtkLmsk);
	field::Field<FLDMASK>  mask  = VTK::createFieldfromVTK<VTKMASK,FLDMASK>(vtkMask);

	if (vecf.get_m() != 3) {
		vtkErrorMacro("Field <"<<this->vecFName<<"> is not vectorial! Aborting...");
		return 0;
	}

	// Recover weights
	VTKARRAY *vtke1, *vtke2, *vtke3;
	if (iscelld) {
		vtke1 = VTKARRAY::SafeDownCast( d->GetCellData()->GetArray("e1") );
		vtke2 = VTKARRAY::SafeDownCast( d->GetCellData()->GetArray("e2") );
		vtke3 = VTKARRAY::SafeDownCast( d->GetCellData()->GetArray("e3") );
	} else {
		vtke1 = VTKARRAY::SafeDownCast( d->GetPointData()->GetArray("e1") );
		vtke2 = VTKARRAY::SafeDownCast( d->GetPointData()->GetArray("e2") );
		vtke3 = VTKARRAY::SafeDownCast( d->GetPointData()->GetArray("e3") );		
	}

	if (vtke1 == NULL || vtke2 == NULL || vtke3 == NULL) {
		vtkErrorMacro("Mesh weights (e1, e2 and e3) need to be loaded to proceed!");
		return 0;
	}

	// Convert to field arrays
	field::Field<FLDARRAY> e1, e2, e3;
	if (vtke1) e1 = VTK::createFieldfromVTK<VTKARRAY,FLDARRAY>(vtke1);
	if (vtke2) e2 = VTK::createFieldfromVTK<VTKARRAY,FLDARRAY>(vtke2);
	if (vtke3) e3 = VTK::createFieldfromVTK<VTKARRAY,FLDARRAY>(vtke3);

	// Create the normals, here they are a vector on the z direction 
	// as of e1t*e2t*[0,0,1]
	field::Field<FLDARRAY> norm(xyz.len(),3,0.);
	if (this->computeCirc) {
		#pragma omp parallel
		{
		for(int ii=OMP_THREAD_NUM; ii<norm.get_n(); ii+=OMP_NUM_THREADS) {
			norm[ii][2] = e1[ii][0]*e2[ii][0];
		}
		}
	}

	// Loop the vortices
	// For each vortex compute the properties and store them
	// in the table.
	PROJ::Projection p;
	int nvortices = (int)(vortexList.size());

	// Create a field for each column of the table and set its value
	// per vortex
	field::Field<FLDARRAY> bcent(nvortices,3), pcent(nvortices,3), size(nvortices,1), rot(nvortices,3), 
		s_abs(nvortices,1), s_rel(nvortices,1), gamma(nvortices,1);

	for (int ivort=0; ivort<nvortices; ++ivort) {
		// Baricenter
		if (this->computeBCent) {
			v3::V3 bxyz = vortexList[ivort].baricenter(xyz);
			bcent[ivort][0] = bxyz[0];
			bcent[ivort][1] = bxyz[1];
			bcent[ivort][2] = bxyz[2];
			// Transform the center to Lon, Lat
			p.transform_point(this->projName,"degrees",bcent[ivort][0],bcent[ivort][1]);	
		}
		// Computation of the center
		if (this->computeCent) {
			int ic = vortexList[ivort].icenter<FLDARRAY>(scaf);
			pcent[ivort][0] = xyz[ic][0];
			pcent[ivort][1] = xyz[ic][1];
			pcent[ivort][2] = xyz[ic][2];
			// Transform the center to Lon, Lat
			p.transform_point(this->projName,"degrees",pcent[ivort][0],pcent[ivort][1]);
		}
		// Computation of the vortex size
		if (this->computeSize) {
			size[ivort][0] = vortexList[ivort].size<FLDARRAY>(xyz,scaf);
		}
		// Rotation axis
		if (this->computeRAxis) {
			rot[ivort][0] = vortexList[ivort].rotation<FLDARRAY>(vecf)[0];
			rot[ivort][1] = vortexList[ivort].rotation<FLDARRAY>(vecf)[1];
			rot[ivort][2] = vortexList[ivort].rotation<FLDARRAY>(vecf)[2];
			// Change the value of the mask to negative to those vortices that are
			// anticyclonic
			// 1: cyclonic, -1: anticyclonic
			for (int ii=0; ii<vortexList[ivort].get_n() && rot[ivort][2] < 0 && this->changemask; ++ii)
				lmask[vortexList[ivort][ii]][0] *= -1.;
		}
		// Vortex absolute strength
		if (this->computeAStr) {
			s_abs[ivort][0] = vortexList[ivort].absolute_strength<FLDARRAY>(vecf);
		}
		// Vortex absolute strength
		if (this->computeRStr) {
			s_rel[ivort][0] = vortexList[ivort].relative_strength<FLDARRAY>(scaf);
		}
		// Circulation
		if (this->computeCirc) {
			gamma[ivort][0] = vortexList[ivort].circulation<FLDARRAY>(norm,vecf);
		}
	}

	// Set the columns to the table
	// Baricenter
	if (this->computeBCent) {
		vtkColumn = VTK::createVTKfromField<VTKARRAY,FLDARRAY>("Baricenter",bcent);
		t->AddColumn(vtkColumn); vtkColumn->Delete();		
	}
	// Vortex center
	if (this->computeCent) {
		vtkColumn = VTK::createVTKfromField<VTKARRAY,FLDARRAY>("Center",pcent);
		t->AddColumn(vtkColumn); vtkColumn->Delete();
	}
	// Size
	if (this->computeSize) {
		vtkColumn = VTK::createVTKfromField<VTKARRAY,FLDARRAY>("Size",size);
		t->AddColumn(vtkColumn); vtkColumn->Delete();
	}
	// Rotation
	if (this->computeRAxis) {
		vtkColumn = VTK::createVTKfromField<VTKARRAY,FLDARRAY>("Rotation",rot);
		t->AddColumn(vtkColumn); vtkColumn->Delete();
	}
	// Absolute strength
	if (this->computeAStr) {
		vtkColumn = VTK::createVTKfromField<VTKARRAY,FLDARRAY>("Abs Str.",s_abs);
		t->AddColumn(vtkColumn); vtkColumn->Delete();
	}
	// Relative strength
	if (this->computeRStr) {
		vtkColumn = VTK::createVTKfromField<VTKARRAY,FLDARRAY>("Rel Str.",s_rel);
		t->AddColumn(vtkColumn); vtkColumn->Delete();
	}
	// Circulation
	if (this->computeCirc) {
		vtkColumn = VTK::createVTKfromField<VTKARRAY,FLDARRAY>("Circulation",gamma);
		t->AddColumn(vtkColumn); vtkColumn->Delete();
	}

	// Set mask array
	vtkLmsk = VTK::createVTKfromField<VTKARRAY,FLDARRAY>("vortex mask",lmask);

	if (iscelld) {
		d->GetCellData()->RemoveArray("vortex mask");
		d->GetCellData()->AddArray(vtkLmsk);
	} else {
		d->GetPointData()->RemoveArray("vortex mask");
		d->GetPointData()->AddArray(vtkLmsk);
	}

	vtkLmsk->Delete();

	return 1;
}

//----------------------------------------------------------------------------
vtkStringArray *vtkOGSVortexDetection::GetProjections() {

	this->Projections->Delete();
	this->Projections = VTK::createVTKstrf("Projections",10,NULL);

	this->Projections->SetValue(0,"Mercator");
	this->Projections->SetValue(1,"Cylindrical");
	this->Projections->SetValue(2,"Google");
	this->Projections->SetValue(3,"Mollweide");
	this->Projections->SetValue(4,"Orthographic");
	this->Projections->SetValue(5,"Robinson");
	this->Projections->SetValue(6,"Satellite");
	this->Projections->SetValue(7,"Eckert IV");
	this->Projections->SetValue(8,"Equal Earth");
	this->Projections->SetValue(9,"EPSG 3857");

	return this->Projections;
}

void vtkOGSVortexDetection::SetProjection(const char *proj) {
	this->projName = std::string(proj);
	this->Modified();
}