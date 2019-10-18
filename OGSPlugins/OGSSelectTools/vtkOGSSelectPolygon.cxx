/*=========================================================================

  Program:   OGSSelectTools
  Module:    vtkOGSSelectPolygon.cxx

  Copyright (c) 2019 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkOGSSelectPolygon.h"

#include "vtkTypeUInt8Array.h"
#include "vtkFloatArray.h"
#include "vtkDoubleArray.h"
#include "vtkCellData.h"
#include "vtkPointData.h"
#include "vtkDataArraySelection.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkUnstructuredGrid.h"

#include "vtkObjectFactory.h"

#include <cstdint>
#include <algorithm>
#include <vector>

#ifdef __linux__
// Include OpenMP when working with GCC
#include <omp.h>
#define OMP_NUM_THREADS omp_get_num_threads()
#define OMP_THREAD_NUM  omp_get_thread_num()
#else
#define OMP_NUM_THREADS 1
#define OMP_THREAD_NUM  0
#endif

#ifdef PARAVIEW_USE_MPI
#include "vtkMultiProcessController.h"
vtkCxxSetObjectMacro(vtkOGSSelectPolygon, Controller, vtkMultiProcessController);
#endif

vtkStandardNewMacro(vtkOGSSelectPolygon);

//----------------------------------------------------------------------------

/*
	Macro to set the array precision 
*/
#define FLDMASK uint8_t
#define VTKMASK vtkTypeUInt8Array

#include "../_utils/V3.h"
#include "../_utils/vtkOperations.hpp"
#include "../_utils/field.h"
#include "../_utils/vtkFields.hpp"
#include "../_utils/geometry.h"
#include "../_utils/Projection.h"

//----------------------------------------------------------------------------
void strsplit(const std::string& str, std::vector<std::string> &cont, char delim) {
    std::size_t current, previous = 0;
    current = str.find(delim);
    while (current != std::string::npos) {
        cont.push_back(str.substr(previous, current - previous));
        previous = current + 1;
        current = str.find(delim, previous);
    }
    cont.push_back(str.substr(previous, current - previous));
}

//----------------------------------------------------------------------------
vtkOGSSelectPolygon::vtkOGSSelectPolygon() {

	this->strPoints  = NULL;
	this->nProcs     = 0;
	this->Invert     = 0;
	this->dfact      = 1000.;

	#ifdef PARAVIEW_USE_MPI
		this->Controller = NULL;
		this->SetController(vtkMultiProcessController::GetGlobalController());
	#endif
}

//----------------------------------------------------------------------------
vtkOGSSelectPolygon::~vtkOGSSelectPolygon() {
	this->SetstrPoints(NULL);

	#ifdef PARAVIEW_USE_MPI
		this->SetController(NULL);	
	#endif
}

//----------------------------------------------------------------------------
int vtkOGSSelectPolygon::RequestData(vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector, vtkInformationVector *outputVector) {
	// Get the info objects
	vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
	vtkInformation *outInfo = outputVector->GetInformationObject(0);

	// Stop all threads except from the master to execute
	#ifdef PARAVIEW_USE_MPI
	if (this->procId > 0) return 1;
	#endif

	// Get the input and output
	vtkDataSet *input = vtkDataSet::SafeDownCast(
		inInfo->Get(vtkDataObject::DATA_OBJECT()));
	vtkUnstructuredGrid *output = vtkUnstructuredGrid::SafeDownCast(
		outInfo->Get(vtkDataObject::DATA_OBJECT()));

	this->UpdateProgress(0.0);

	// Obtain information on the projection (Metadata array)
	vtkStringArray *vtkmetadata = vtkStringArray::SafeDownCast(
		input->GetFieldData()->GetAbstractArray("Metadata"));
	this->dfact    = (vtkmetadata != NULL) ? std::stod( vtkmetadata->GetValue(2) ) : this->dfact;
	this->projName = (vtkmetadata != NULL) ? vtkmetadata->GetValue(7) : std::string("Mercator");

	// Parse the points inputted by the user, project to meters and define the polygon
	std::vector<std::string> aux;
	strsplit(std::string(this->strPoints),aux,'\n'); // Split by endline

	std::vector<Geom::Point<double>> points;
	PROJ::Projection p;	

	// Loop the user inputed points
	for (std::string str : aux) {
		// Split again by ;
		std::vector<std::string> aux2;
		strsplit(str,aux2,' ');
		// Check correctness
		if (str.size() == 0 || aux2.size() != 2) {
			vtkErrorMacro("Problems parsing points! Check the input points.");
			return 0;
		}
		// Convert to double
		double lon = std::stod(aux2[1]), lat = std::stod(aux2[0]);
		// Project
		std::transform(this->projName.begin(), this->projName.end(), this->projName.begin(), ::tolower);
		p.transform_point("degrees",this->projName,lon,lat);
		// Store point
		points.push_back( Geom::Point<double>(lon,lat,0.) );
	}

	// Define the polygon and compute the bounding box
	Geom::Polygon<double> poly((int)(points.size()),points.data());
	
	// Compute the polygon bounding box for a faster performance
	Geom::Ball<double> bbox(poly);
	poly.set_bbox(bbox);

	this->UpdateProgress(0.1);

	// Understand whether we are under cell or point data and compute the points of
	// the mesh
	int n_cell_vars  = input->GetCellData()->GetNumberOfArrays();
	int n_point_vars = input->GetPointData()->GetNumberOfArrays();

	bool iscelld = (n_cell_vars > n_point_vars) ? true : false;

	v3::V3v xyz = (iscelld) ? VTK::getVTKCellCenters(input,this->dfact) : VTK::getVTKCellPoints(input,this->dfact);

	// Generate a new field (initialized at zero) that will be used as cutting mask
	field::Field<FLDMASK> cutmask(xyz.len(),1);

	this->UpdateProgress(0.2);

	// Loop and update cutting mask (Mesh loop, can be parallelized)
	#pragma omp parallel shared(cutmask,poly,xyz)
	{
	for (int ii = 0 + OMP_THREAD_NUM; ii < cutmask.get_n(); ii += OMP_NUM_THREADS) {
		// Check if the point of the mesh is inside the polygon and set the cutmask
		// accordingly
		cutmask[ii][0] = (poly > Geom::Point<double>(xyz[ii][0],xyz[ii][1],0.)) ? 1 : 0;
	}
	}

	// Convert field to vtkArray and add it to input
	VTKMASK *vtkcutmask;
	vtkcutmask = VTK::createVTKfromField<VTKMASK,FLDMASK>("CutMask",cutmask);

	if (iscelld) {
		input->GetCellData()->AddArray(vtkcutmask);
		// Force to use the CutMask to produce the Threshold
		this->Superclass::SetInputArrayToProcess(0,0,0,
			vtkDataObject::FIELD_ASSOCIATION_CELLS,"CutMask");
	} else {
		input->GetPointData()->AddArray(vtkcutmask);
		// Force to use the CutMask to produce the Threshold
		this->Superclass::SetInputArrayToProcess(0,0,0,
			vtkDataObject::FIELD_ASSOCIATION_POINTS,"CutMask");
	}

	this->UpdateProgress(0.4);

	// Force ThresholdBetween to obtain values that are greater than 0
	if (this->Invert)
		this->Superclass::ThresholdBetween(0.,0.5);
	else
		this->Superclass::ThresholdBetween(0.5,1.);

	this->UpdateProgress(0.6);

	// Run the actual threshold filter
	this->Superclass::RequestData(NULL,inputVector,outputVector);

	this->UpdateProgress(0.8);

	// Cleanup the output by deleting the CutMask and the basins mask	
	if (iscelld)
		output->GetCellData()->RemoveArray("CutMask");
	else 
		output->GetPointData()->RemoveArray("CutMask");

	vtkcutmask->Delete();


	// Return
	this->UpdateProgress(1.0);
	return 1;
}