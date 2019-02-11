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
#include "../_utils/OGS.hpp"

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
	// Get the name of the master file to read
	vtkSetStringMacro(FileName);
	vtkGetStringMacro(FileName);

	// Description:
	// If false, do not include the meshmask. True by default.
	vtkGetMacro(RMeshMask, int);
	vtkSetMacro(RMeshMask, int);
	vtkBooleanMacro(RMeshMask, int);

	// Description:
	// Lets the user select a multiplier factor for the depth
	vtkGetMacro(DepthScale, double);
	vtkSetMacro(DepthScale, double);

	// Description:
	// The following methods allow selective reading of the mask variables.
	// By default, ALL variables are read, but this can be modified 
	// (e.g. from the ParaView GUI).
	int GetNumberOfMaskArrays();
	const char * GetMaskArrayName(int index);
	int GetMaskArrayIndex(const char* name);
	int GetMaskArrayStatus(const char *name);
	void SetMaskArrayStatus(const char* name, int status);
	void DisableAllMaskArrays();
	void EnableAllMaskArrays();

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

	// Description:
	// The following methods allow selective reading of the forcing variables.
	// By default, ALL variables are read, but this can be modified 
	// (e.g. from the ParaView GUI).
	int GetNumberOfForcingArrays();
	const char * GetForcingArrayName(int index);
	int GetForcingArrayIndex(const char* name);
	int GetForcingArrayStatus(const char *name);
	void SetForcingArrayStatus(const char* name, int status);
	void DisableAllForcingArrays();
	void EnableAllForcingArrays();

	// Description:
	// The following methods allow selective reading of the general variables.
	// By default, ALL variables are read, but this can be modified 
	// (e.g. from the ParaView GUI).
	int GetNumberOfGeneralArrays();
	const char * GetGeneralArrayName(int index);
	int GetGeneralArrayIndex(const char* name);
	int GetGeneralArrayStatus(const char *name);
	void SetGeneralArrayStatus(const char* name, int status);
	void DisableAllGeneralArrays();
	void EnableAllGeneralArrays();

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

	int RMeshMask;
	double DepthScale;

	vtkDataArraySelection* MaskDataArraySelection;
	vtkDataArraySelection* AvePhysDataArraySelection;
	vtkDataArraySelection* AveFreqDataArraySelection;
	vtkDataArraySelection* ForcingDataArraySelection;
	vtkDataArraySelection* GeneralDataArraySelection;

	#ifdef PARAVIEW_USE_MPI
		vtkMultiProcessController* Controller;
	#endif

private:
	vtkOGSReader(const vtkOGSReader&) = delete;
	void operator=(const vtkOGSReader&) = delete;

	vtkRectilinearGrid* Mesh;

	ogs::OGS ogsdata;

	int abort;

	void DeleteMesh();
};

#endif