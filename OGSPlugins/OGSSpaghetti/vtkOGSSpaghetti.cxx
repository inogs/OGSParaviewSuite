/*=========================================================================

  Program:   OGSSpaghetti
  Module:    vtkOGSSpaghetti.cxx

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkOGSSpaghetti.h"

#include "vtkCell.h"
#include "vtkCellData.h"
#include "vtkCellLocator.h"
#include "vtkStringArray.h"
#include "vtkGenericCell.h"
#include "vtkPointData.h"
#include "vtkPointSet.h"
#include "vtkDataSet.h"
#include "vtkTable.h"
#include "vtkDataArray.h"
#include "vtkTypeUInt8Array.h"
#include "vtkFloatArray.h"
#include "vtkDoubleArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkStaticCellLocator.h"
#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#include "vtkObjectFactory.h"

#include <cstdint>
#include <cmath>
#include <ctime>
#include <chrono>
#include <string>

vtkStandardNewMacro(vtkOGSSpaghetti);
vtkCxxSetObjectMacro(vtkOGSSpaghetti, CellLocatorPrototype, vtkAbstractCellLocator);

#ifdef PARAVIEW_USE_MPI
#include "vtkMultiProcessController.h"
vtkCxxSetObjectMacro(vtkOGSSpaghetti, Controller, vtkMultiProcessController);
#endif

#define CELL_TOLERANCE_FACTOR_SQR 1e-6

//----------------------------------------------------------------------------

/*
	Macro to set the array precision 
*/
#define FLDARRAY double
#define VTKARRAY vtkDoubleArray
#define FLDMASK uint8_t
#define VTKMASK vtkTypeUInt8Array

#include "../_utils/vtkFields.hpp"
#include "../_utils/OGS.hpp"
#include "../_utils/netcdfio.hpp"
#include "../_utils/fieldOperations.hpp"
#include "../_utils/vtkOperations.hpp"

void strsplit(std::string str, std::string splitBy, std::vector<std::string> &tokens) {
    // Store the original string in the array, so we can loop the rest of the algorithm.
    tokens.push_back(str);

    // Store the split index in a 'size_t' (unsigned integer) type.
    size_t splitAt;
    // Store the size of what we're splicing out.
    size_t splitLen = splitBy.size();
    // Create a string for temporarily storing the fragment we're processing.
    std::string frag;
    // Loop infinitely - break is internal.
    while(true) {
        // Store the last string in the vector, which is the only logical candidate for processing.
        frag = tokens.back();
        // The index where the split is.
        splitAt = frag.find(splitBy);
        // If we didn't find a new split point...
        if(splitAt == std::string::npos)
            break; // Break the loop and (implicitly) return.
        // Put everything from the left side of the split where the string being processed used to be.
        tokens.back() = frag.substr(0, splitAt);
        // Push everything from the right side of the split to the next empty index in the vector.
        tokens.push_back(frag.substr(splitAt+splitLen, frag.size()-(splitAt+splitLen)));
    }
}

//----------------------------------------------------------------------------
void RecoverMasterFileName(std::string &fname, vtkDataSet *input) {

	// Recover the master file name from the metadata array
	// Return whether we need to stop executing or not

	vtkStringArray *vtkmetadata = vtkStringArray::SafeDownCast(
		input->GetFieldData()->GetAbstractArray("Metadata"));

	// If successful, recover the file name
	if (vtkmetadata)
		fname = vtkmetadata->GetValue(4);
}

//----------------------------------------------------------------------------
int ComputeCellId(vtkDataSet *input, vtkDataSet *source) {
	int cId = -1;

	// We need the weights as they are part of a function
	double *weights; weights = new double[source->GetMaxCellSize()];

	// Create the cell locator object
	vtkCellLocator *cellLocator = vtkCellLocator::New();
	cellLocator->SetDataSet(source); cellLocator->BuildLocator();

	// Get the xyz coordinate of the point in the input dataset
	// then, find the cell id that contains xyz
	vtkIdType cellId = 0; int subId = 0;
	double dist2;
	v3::V3 xyz, pcoords, closestPoint;
	vtkNew<vtkGenericCell> gcell;
	
	input->GetPoint(0, &xyz[0]);
	cellLocator->FindClosestPoint(&xyz[0], &closestPoint[0], gcell.GetPointer(), cellId, subId, dist2);

	// Evaluate interpolation weights
	if (cellId >= 0) {
		// Compute a tolerance proportional to the cell length.
		gcell->EvaluatePosition(&xyz[0], &closestPoint[0], subId, &pcoords[0], dist2, weights);
		// Abort if the distance is too big
		if (dist2 > (gcell->GetLength2() * CELL_TOLERANCE_FACTOR_SQR))
			return cId;
		// Store the cell id
		cId = (int)(cellId);
	}

	return cId;
}

//----------------------------------------------------------------------------
void ComputeCellIds(std::vector<int> &cellIds, vtkDataSet *input, vtkDataSet *source) {

	// We need the weights as they are part of a function
	double *weights; weights = new double[source->GetMaxCellSize()];

	// Create the cell locator object
	vtkCellLocator *cellLocator = vtkCellLocator::New();
	cellLocator->SetDataSet(source); cellLocator->BuildLocator();

	// Recover the number of points in the input
	int npoints = input->GetNumberOfPoints();
	cellIds.resize(npoints,-1);

	// Loop the number of points
	vtkNew<vtkGenericCell> gcell;
	for (int pId = 0; pId < npoints; ++pId) {
		// Get the xyz coordinate of the point in the input dataset
		// then, find the cell id that contains xyz
		vtkIdType cellId = 0; int subId = 0;
		double dist2;
		v3::V3 xyz, pcoords, closestPoint; 
		
		input->GetPoint(pId, &xyz[0]);
		cellLocator->FindClosestPoint(&xyz[0], &closestPoint[0], gcell.GetPointer(), cellId, subId, dist2);

		// Evaluate interpolation weights
		if (cellId >= 0) {
			// Compute a tolerance proportional to the cell length.
			gcell->EvaluatePosition(&xyz[0], &closestPoint[0], subId, &pcoords[0], dist2, weights);
			// Abort if the distance is too big
			if (dist2 > (gcell->GetLength2() * CELL_TOLERANCE_FACTOR_SQR))
				continue;
			// Store the cell id
			cellIds[pId] = (int)(cellId);
		}
	}
	delete [] weights;
	cellLocator->Delete();
}

//----------------------------------------------------------------------------
vtkOGSSpaghetti::vtkOGSSpaghetti() {
	this->SetNumberOfInputPorts(2);
	this->SetNumberOfOutputPorts(1);

	this->PointList   = nullptr;
	this->CellList    = nullptr;
	this->field       = nullptr;
	this->field       = NULL;
	this->FolderName  = NULL;
	this->bmask_field = NULL;
	this->cmask_field = NULL;
	this->ii_start    = 0;
	this->ii_end      = 0;
	this->average     = 0;
	this->sId         = 0;
	this->per_coast   = 0;
	this->procId      = 0;
	this->nProcs      = 0;
	this->isReqInfo   = false;
	this->epsi        = 1.e-3;
	this->dfact       = 1000.;

	this->CellLocatorPrototype = nullptr;
	this->TimeValues = vtkStringArray::New();

	#ifdef PARAVIEW_USE_MPI
	this->Controller = nullptr;
	this->SetController(vtkMultiProcessController::GetGlobalController());
	#endif
}

//----------------------------------------------------------------------------
vtkOGSSpaghetti::~vtkOGSSpaghetti() {
	delete this->PointList;
	delete this->CellList;

	this->Setfield(NULL);
	this->SetFolderName(NULL);
	this->Setbmask_field(NULL);
	this->Setcmask_field(NULL);

	this->vtkOGSSpaghetti::SetCellLocatorPrototype(nullptr);
	this->TimeValues->Delete();

	#ifdef PARAVIEW_USE_MPI
	this->SetController(nullptr);
	#endif
}

int vtkOGSSpaghetti::FillInputPortInformation(int vtkNotUsed(port), vtkInformation *info) {
  info->Remove(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE());
  info->Append(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
  return 1;
}

//----------------------------------------------------------------------------
int vtkOGSSpaghetti::RequestInformation(vtkInformation* vtkNotUsed(request),
  vtkInformationVector** inputVector, vtkInformationVector* outputVector) {

  	// Grab the input
	vtkInformation *outInfo = outputVector->GetInformationObject(0);
	vtkInformation *srcInfo = inputVector[1]->GetInformationObject(0);

	// Handle the parallel controller
	#ifdef PARAVIEW_USE_MPI
	if (this->Controller->GetNumberOfProcesses() > 1) {
		this->procId = this->Controller->GetLocalProcessId();
		this->nProcs = this->Controller->GetNumberOfProcesses();
	}
	#endif

	// Populate the TimeValues array. For that we will use the data
	// stored in the timesteps to be consistent and since the metadata
	// array might not be available from the beginning.
	int ntsteps = srcInfo->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
	double *tsteps; tsteps = srcInfo->Get(vtkStreamingDemandDrivenPipeline::TIME_STEPS());

	char buff[256];
	this->TimeValues->Delete();
	this->TimeValues = VTK::createVTKstrf("TimeValues",ntsteps,NULL);
	
	for (int ii = 0; ii < ntsteps; ++ii) {
		// Convert to struct tm
		time_t time = std::chrono::system_clock::to_time_t(std::chrono::system_clock::time_point(
			std::chrono::duration_cast<std::chrono::seconds>(std::chrono::duration<double>(tsteps[ii]))));
		struct tm tm = *localtime(&time);

		// Format and display
		strftime(buff,256,"%Y%m%d-%H:%M:%S",&tm);
		this->TimeValues->SetValue(ii,buff);
	}

	// The output data of this filter has no time associated with it. It is the
	// result of computations that happen over all time.
	outInfo->Remove(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
	outInfo->Remove(vtkStreamingDemandDrivenPipeline::TIME_RANGE());

	return 1;
}

//----------------------------------------------------------------------------
int vtkOGSSpaghetti::RequestData(vtkInformation *request, 
	vtkInformationVector **inputVector, vtkInformationVector *outputVector) {
	
	// get the info objects
	vtkInformation *inInfo  = inputVector[0]->GetInformationObject(0);
	vtkInformation *srcInfo = inputVector[1]->GetInformationObject(0);
	vtkInformation *outInfo = outputVector->GetInformationObject(0);

	int ntsteps = srcInfo->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS());

	// input contains the interpolating line information (number of points, etc)
	vtkDataSet *input = vtkDataSet::SafeDownCast(
		inInfo->Get(vtkDataObject::DATA_OBJECT()));
	// Source contains the mesh where to interpolate from
	vtkDataSet *source = vtkDataSet::SafeDownCast(
		srcInfo->Get(vtkDataObject::DATA_OBJECT()));
	// Output is a vtkTable with the interpolated data per each timestep
	vtkTable *output = vtkTable::SafeDownCast(
		outInfo->Get(vtkDataObject::DATA_OBJECT()));

	// Check input field
	if (std::string(this->field) == std::string("coast mask"))  { vtkErrorMacro("Wrong input field! Aborting..."); return 0; }
	if (std::string(this->field) == std::string("basins mask")) { vtkErrorMacro("Wrong input field! Aborting..."); return 0; }
	if (std::string(this->field) == std::string("e1"))          { vtkErrorMacro("Wrong input field! Aborting..."); return 0; }
	if (std::string(this->field) == std::string("e2"))          { vtkErrorMacro("Wrong input field! Aborting..."); return 0; }
	if (std::string(this->field) == std::string("e3"))          { vtkErrorMacro("Wrong input field! Aborting..."); return 0; }
	if (std::string(this->field) == std::string("e1"))          { vtkErrorMacro("Wrong input field! Aborting..."); return 0; }

	// Decide which algorithm to use
	int retval = 0;
	if (this->average)
		return this->SpaghettiAverage(ntsteps,input,source,output);
	else
		return this->Spaghetti3DDataset(input,source,output);
}

//----------------------------------------------------------------------------
int vtkOGSSpaghetti::Spaghetti3DDataset(vtkDataSet *input, vtkDataSet *source, vtkTable *output) {
	/* INITIALIZATION PHASE

		Master computes the cell list corresponding to the Hovmoller interpolating line.
		Recover (and broadcast) the file name of the OGS master file.

		Initialize arrays.

	*/
	std::string FileName;
	int cellList, projId = 0;
	vtkStringArray *vtkmetadata;

	#ifdef PARAVIEW_USE_MPI

	// Thread 0 has all the information in input and output, therefore is the
	// only computing the array which will later be broadcasted to all ranks
	if (this->procId == 0) {
		// Recover metadata array
		vtkmetadata = vtkStringArray::SafeDownCast(source->GetFieldData()->GetAbstractArray("Metadata"));
		// Recover projection Id
		projId = std::stod( vtkmetadata->GetValue(7) );
		// Compute the cell list corresponding to the Hovmoeller interpolating line.
		cellList = ComputeCellId(input,source);
		// Recover the master file name from source
		RecoverMasterFileName(FileName, source);
	}

	// Broadcast the information to all other threads if applicable
	if (this->nProcs > 1) {
		// Broadcast projection ID
		this->Controller->Broadcast(&projId,1,0);
		// Broadcast master file name
		int str_len = FileName.length();
		this->Controller->Broadcast(&str_len,1,0);
		char buff[128] = ""; sprintf(buff,"%s",FileName.c_str());
		this->Controller->Broadcast(buff,str_len,0);
		FileName = std::string(buff);
		// Broadcast cellList
		this->Controller->Broadcast(&cellList,1,0);
	}

	#else

	// Recover metadata array
	vtkmetadata = vtkStringArray::SafeDownCast(source->GetFieldData()->GetAbstractArray("Metadata"));
	// Recover projection Id
	projId = std::stod( vtkmetadata->GetValue(7) );

	// This is the normal non-parallel algorithm
	cellList = ComputeCellId(input,source);
	RecoverMasterFileName(FileName, source);

	#endif

	// Read the OGS file to be able to generate the paths to the variable files
	ogs::OGS ogsdata(FileName);
	if (ogsdata.readMainFile() == -1) { 
		vtkErrorMacro("Cannot read <"<<FileName<<">!\nAborting");
		return 0;
	}

	// Read the mesh data (necessary to load the files)
	if (ogsdata.readMesh(projId) < 0) {
		vtkErrorMacro("Problems reading the mesh!\nAborting.");
		return 0;	
	}

	this->UpdateProgress(0.1);

	/* TRAVERSE PHASE

		Open and read the NetCDF files containing the variables that we previously
		defined. Store the data in a master field for posterior reduction.

	*/
	int time_interval[2] = {this->ii_start,this->ii_end};

	// Parallel partition
	#ifdef PARAVIEW_USE_MPI

	if (this->nProcs > 1) {
		// Main thread (0) contains values for ii_start and ii_end
		// They must be sent to the other processes
		this->Controller->Broadcast(time_interval,2,0);

		// Now everyone should have the range of start and end times
		// We must split equally among the threads. We must also control
		// that the number of threads is less or equal than the intervals
		// requested.
		int range = time_interval[1] - time_interval[0];
		if (this->nProcs < range) {
			// We split normally among processes assuming no remainder
			int rangePerProcess = std::floor(range/this->nProcs);
			this->ii_start = time_interval[0] + this->procId*rangePerProcess;
			this->ii_end   = this->ii_start + rangePerProcess;
			// Handle the remainder
			int remainder = range - rangePerProcess*this->nProcs;
			if (remainder > this->procId){
				this->ii_start += this->procId;
				this->ii_end   += this->procId + 1;
			} else {
				this->ii_start += remainder;
				this->ii_end   += remainder;
			}
		} else {
			// Each process will forcefully conduct one instant.
			this->ii_start = (this->procId < time_interval[1]) ? this->procId     : time_interval[1];
			this->ii_end   = (this->procId < time_interval[1]) ? this->procId + 1 : time_interval[1];
		}
	}

	#endif

	// Ensure the validity of the time range
	if (time_interval[0] > time_interval[1]) {
		if (this->procId == 0) 
			vtkErrorMacro("End time is greater than initial time! Please select a valid time range.");
		return 0;
	}

	// Loop the instants, prepare arrays and variables for the computation of the mean
	field::Field<FLDARRAY> arrayTemp, spag_array(time_interval[1]-time_interval[0],1,0.);
	
	// For each timestep
	for(int ii = this->ii_start; ii < this->ii_end; ++ii) {
		// Load the variable selected by the user
		if (std::string(this->field) == std::string("Velocity")) {
			vtkErrorMacro("Cannot perform the Hovmoeller plot of the Velocity for now..."); return 0;
		}

		// Load the variable on a temporal field
		arrayTemp.set_dim(ogsdata.ncells(),1);

		if ( NetCDF::readNetCDF(ogsdata.var_path(this->field,ii).c_str(), 
			ogsdata.var_vname(this->field), arrayTemp) != NETCDF_OK ) {
			vtkErrorMacro("Cannot read variable <"<<this->field<<"> in NetCDF! Aborting!"); return 0;
		}

		// Fill the Spaghetti array using the precomputed cell list
		spag_array[ii-time_interval[0]][0] = (cellList > 0) ? arrayTemp[cellList][0] : 0.;

		arrayTemp.clear();
		this->UpdateProgress(0.1+0.7/(this->ii_end-this->ii_start)*(ii - this->ii_start));
	}

	/* REDUCTION PHASE

		If run in parallel with more than one rank, 
		the reduction of the fields is done here.

	*/
	#ifdef PARAVIEW_USE_MPI

	// Only reduce if we have more than 1 process
	if (this->nProcs > 1) {

		field::Field<FLDARRAY> spag_array_tmp(time_interval[1]-time_interval[0],1,0.);

		// Reduce the Spaghetti array
		this->Controller->Reduce(spag_array.data(),spag_array_tmp.data(),spag_array.get_sz(),
			vtkCommunicator::StandardOperations::SUM_OP,0);

		spag_array.set_val( spag_array_tmp.data() ); spag_array_tmp.clear();
	}

	#endif

	this->UpdateProgress(0.9);


	/* FINALIZE

		Finalization, the master process stores the arrays
		inside the output.

	*/
	// Build output table
	if (this->procId == 0) {
		// Recover metadata array
/*		vtkStringArray *vtkmetadata = vtkStringArray::SafeDownCast(
			source->GetFieldData()->GetAbstractArray("Metadata"));*/
		// Recover datevec
		std::string datevec = vtkmetadata->GetValue(1);
		std::vector<std::string> vdatevec;
		strsplit(datevec,";",vdatevec);

		// Auxiliar field for columns
		int npoints = input->GetNumberOfPoints();
		field::Field<FLDARRAY> column(npoints,1);
		VTKARRAY *vtkColumn;

		// Write the time instants
		for (int ii = time_interval[0]; ii < time_interval[1]; ++ii) {
			column.set_val( spag_array[ii-time_interval[0]] );
			vtkColumn = VTK::createVTKfromField<VTKARRAY,FLDARRAY>(vdatevec[ii+1],column);
			output->AddColumn(vtkColumn); vtkColumn->Delete();
		}
	}

	this->UpdateProgress(1.);
	return 1;
}

//----------------------------------------------------------------------------
int vtkOGSSpaghetti::SpaghettiAverage(int ntsteps, vtkDataSet *input, vtkDataSet *source, vtkTable *output) {

	// There is no need to parallelize this algorithm
	#ifdef PARAVIEW_USE_MPI
	if (this->procId > 0) return 1;
	#endif

	// Recover Metadata array (depth factor)
	vtkStringArray *vtkmetadata = vtkStringArray::SafeDownCast(
		source->GetFieldData()->GetAbstractArray("Metadata"));

	// This section is only executed once, to populate the xyz and
	// cId2zId arrays. Successive iterations should not execute.
	// This section is included here since RequestInformation gives
	// troubles when restarting.
	if (this->xyz.isempty() || this->isReqInfo) {

		this->isReqInfo = false;
		
		this->dfact = (vtkmetadata != NULL) ? std::stod( vtkmetadata->GetValue(2) ) : 1000.;
		
		if (vtkmetadata == NULL) 
			vtkWarningMacro("Field array Metadata not found! Depth factor set to 1000. automatically.");

		// Recover cell or point coordinates
		this->xyz = VTK::getVTKCellCenters(source,this->dfact);

		// Up to this point we have the cell centers or point coordinates correctly
		// stored under "xyz". Now we shall find the number of unique z coordinates or,
		// depending on the user input, the coordinates of each depth level, as well as
		// its mesh connectivity (cId2zId).
		this->cId2zId = field::countDepthLevels(this->xyz,this->zcoords,this->epsi,false);
	}

	this->UpdateProgress(.1);

	// Recover the basins and coasts mask
	// and add them to the output
	int nbasins = 21, ncoasts = 3, nStat = 9;
	VTKMASK *vtkMask;
	field::Field<FLDMASK> bmask, cmask;
	// Basins mask
	vtkMask = VTKMASK::SafeDownCast( source->GetCellData()->GetArray(this->bmask_field) );
	bmask   = VTK::createFieldfromVTK<VTKMASK,FLDMASK>(vtkMask);
	// Coast mask
	vtkMask = VTKMASK::SafeDownCast( source->GetCellData()->GetArray(this->cmask_field) );
	cmask   = VTK::createFieldfromVTK<VTKMASK,FLDMASK>(vtkMask);

	// Load the mesh weights
	VTKARRAY *vtke3 = NULL;
	vtke3 = VTKARRAY::SafeDownCast( source->GetCellData()->GetArray("e3") );

	if (vtke3 == NULL) {
		vtkErrorMacro("Mesh weights (e1, e2 and e3) need to be loaded to proceed!");
		return 0;
	}

	// Convert to field arrays
	field::Field<FLDARRAY> e3;
	e3 = VTK::createFieldfromVTK<VTKARRAY,FLDARRAY>(vtke3);

	// Load variable stat profile
	field::Field<FLDARRAY> statArray(ntsteps*nbasins*ncoasts*this->zcoords.size()*nStat,1,0.);
	std::string filename = std::string(this->FolderName) + std::string("/") + std::string(this->field) + std::string(".nc");

	if ( NetCDF::readNetCDF(filename.c_str(),this->field,statArray) != NETCDF_OK ) {
		// If file cannot be read or variable does not exist
		vtkErrorMacro("File <"<<filename.c_str()<<"> or variable <"<<this->field<<"> cannot be read!");
		return 0;
	}
	// Retrieve cellIds from interpolating line
	std::vector<int> cellIds;
	ComputeCellIds(cellIds,input,source);

	// Recover datevec
	std::string datevec = vtkmetadata->GetValue(1);
	std::vector<std::string> vdatevec;
	strsplit(datevec,";",vdatevec);

	// Compute the mesh dependent position that does not depend on time
	// we only need to do this once
	std::vector<int> pos1; pos1.resize(cellIds.size(),0);
	double weight_sum = 0.;
	for (int cc = 0; cc < cellIds.size(); ++cc) {
		// Recover cell id
		int cellId = cellIds[cc];
		// Depth index
		int zId = this->cId2zId[cellId][0];
		// In which basin are we? (we need to loop the basins and find which is true)
		int  bId = -1; 
		bool isbasin = false;
		for (bId = 0; bId < bmask.get_m(); ++bId) {
			if (bmask[cellId][bId]) { isbasin = true; break; }
		}
		// In which coast are we?
		int cId = this->per_coast ? cmask[cellId][0] - 1 : 2;
		// Compute the position
		int p = ncoasts*this->zcoords.size()*nStat*bId  + 
		        this->zcoords.size()*nStat*cId          + 
		        nStat*zId                         +
		        this->sId;
		p = (bId > 0 || cId > 0) ? p : -1;
		// Store the position and the sum of weights
		pos1[cc]    = p;
		weight_sum += e3[cellId][0];
	}

	this->UpdateProgress(.25);

	// For each time instant, loop on the cell list and load the data into
	// the table
	field::Field<FLDARRAY> column(1,1);
	VTKARRAY *vtkColumn;
	for (int ii = ii_start; ii < ii_end; ii += 1) {
		column[0][0] = 0.;
		for (int cc = 0; cc < cellIds.size(); ++cc) {
			// Recover ids
			int p      = pos1[cc];
			int cellId = cellIds[cc];
			// Compute position
			int pos = nbasins*ncoasts*this->zcoords.size()*nStat*ii + p;
			// Retrieve value from array
			column[0][0] += (p > 0) ? e3[cellId][0]*statArray[pos][0] : 0.;
		}
		column[0][0] /= weight_sum;
		vtkColumn = VTK::createVTKfromField<VTKARRAY,FLDARRAY>(vdatevec[ii+1],column);
		output->AddColumn(vtkColumn); vtkColumn->Delete();

		this->UpdateProgress(0.25+0.75/(this->ii_end-this->ii_start)*(ii - this->ii_start));
	}

	this->UpdateProgress(1.);
	return 1;
}

//----------------------------------------------------------------------------
void vtkOGSSpaghetti::SetStartTime(const char *tstep) {
	// Obtain the timestep index
	for (int ii = 0; ii < this->TimeValues->GetNumberOfTuples(); ++ii) {
		if (std::string(this->TimeValues->GetValue(ii)) == std::string(tstep)) {
			this->ii_start = ii; break;
		}
	}
	this->Modified();
}

void vtkOGSSpaghetti::SetEndTime(const char *tstep) {
	// Obtain the timestep index
	for (int ii = 0; ii < this->TimeValues->GetNumberOfTuples(); ++ii) {
		if (std::string(this->TimeValues->GetValue(ii)) == std::string(tstep)) {
			this->ii_end = ii; break;
		}
	}
	this->Modified();
}

//----------------------------------------------------------------------------
vtkStringArray *vtkOGSSpaghetti::GetTimeValues() {
	return this->TimeValues;
}

//----------------------------------------------------------------------------
void vtkOGSSpaghetti::SetSourceConnection(vtkAlgorithmOutput* algOutput) {
	this->SetInputConnection(1, algOutput);
}

void vtkOGSSpaghetti::SetSourceData(vtkDataObject *input) {
	this->SetInputData(1, input);
}

vtkDataObject *vtkOGSSpaghetti::GetSource() {
	if (this->GetNumberOfInputConnections(1) < 1)
		return nullptr;

	return this->GetExecutive()->GetInputData(1, 0);
}