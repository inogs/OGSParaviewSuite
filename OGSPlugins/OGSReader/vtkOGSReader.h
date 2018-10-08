// -*- c++ -*-
/*=========================================================================

  Program:   OGSReader
  Module:    vtkOGSReader.h

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef vtkOGSReader_h
#define vtkOGSReader_h

#include "vtkDataSet.h"
#include "vtkRectilinearGridAlgorithm.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#include "vtkPVConfig.h" // For PARAVIEW_USE_MPI

#include "../_utils/OGSdefs.h"

class vtkDataSet;
class vtkDataArraySelection;
class vtkCallbackCommand;

#ifdef PARAVIEW_USE_MPI
	class vtkMultiProcessController;
#endif

//----------------------------------------------------------------------------

// Filter class
class vtkOGSReader : public vtkRectilinearGridAlgorithm {
public:
	static vtkOGSReader* New();
	vtkTypeMacro(vtkOGSReader, vtkRectilinearGridAlgorithm);

	// Description:
	// Get the number of data fields at the cell centers.
	vtkGetMacro(NumberOfAvePhysFields, int);
	vtkGetMacro(NumberOfAveFreqFields, int);

	// Description:
	// Get the number of data components at the nodes and cells.
	vtkGetMacro(NumberOfAvePhysComponents, int);
	vtkGetMacro(NumberOfAveFreqComponents, int);

	// Description:
	// Get the name of the master file to read
	vtkSetStringMacro(FileName);
//	vtkGetStringMacro(FileName);

	// Description:
	// If false, do not include the sub basins. True by default.
	vtkGetMacro(SubBasinsMask, int);
	vtkSetMacro(SubBasinsMask, int);
	vtkBooleanMacro(SubBasinsMask, int);

	// Description:
	// If false, do not include the sub basins. True by default.
	vtkGetMacro(CoastsMask, int);
	vtkSetMacro(CoastsMask, int);
	vtkBooleanMacro(CoastsMask, int);

	// Description:
	// Lets the user select a multiplier factor for the depth
	vtkGetMacro(DepthScale, double);
	vtkSetMacro(DepthScale, double);

	// Description:
	// The following methods allow selective reading of the physical variables.
	// By default, ALL variables are read, but this can be modified 
	// (e.g. from the ParaView GUI).
	int GetNumberOfAvePhysArrays();
	const char * GetAvePhysArrayName(int index);
	int GetAvePhysArrayIndex(const char* name);
	int GetAvePhysArrayStatus(const char *name);
	void SetAvePhysArrayStatus(const char* name, int status);
	void DisableAllAvePhysArrays();
	void EnableAllAvePhysArrays();

	// Description:
	// The following methods allow selective reading of the biogeochemical variables.
	// By default, ALL variables are read, but this can be modified 
	// (e.g. from the ParaView GUI).
	int GetNumberOfAveFreqArrays();
	const char * GetAveFreqArrayName(int index);
	int GetAveFreqArrayIndex(const char* name);
	int GetAveFreqArrayStatus(const char *name);
	void SetAveFreqArrayStatus(const char* name, int status);
	void DisableAllAveFreqArrays();
	void EnableAllAveFreqArrays();

	#ifdef PARAVIEW_USE_MPI
		// Description:
		// Set the controller use in compositing (set to
		// the global controller by default)
		// If not using the default, this must be called before any
		// other methods.
		virtual void SetController(vtkMultiProcessController* controller);
	#endif

protected:
	vtkOGSReader();
	~vtkOGSReader() override;
	int RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector*) VTK_OVERRIDE;
	int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*) VTK_OVERRIDE;

	// Variables
	char *FileName;

	int SubBasinsMask, CoastsMask;
	double DepthScale;

	int NumberOfAvePhysFields, NumberOfAvePhysComponents;
	int NumberOfAveFreqFields, NumberOfAveFreqComponents;

	vtkDataArraySelection* AvePhysDataArraySelection;
	vtkDataArraySelection* AveFreqDataArraySelection;

	#ifdef PARAVIEW_USE_MPI
		vtkMultiProcessController* Controller;
	#endif

private:
	vtkOGSReader(const vtkOGSReader&) = delete;
	void operator=(const vtkOGSReader&) = delete;

	char meshfile[512];

	vtkRectilinearGrid* Mesh;

	ave_var  ave_phys, ave_freq;
	ogs_time timeStepInfo;

	void DeleteMesh();
};

#endif