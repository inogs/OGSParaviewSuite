// -*- c++ -*-
/*=========================================================================

  Program:   OGSSpatialStats
  Module:    vtkOGSSpatialStatsFromFile.h

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef vtkOGSSpatialStatsFromFile_h
#define vtkOGSSpatialStatsFromFile_h

#include "vtkDataSet.h"
#include "vtkRectilinearGridAlgorithm.h"

#include "vtkPVConfig.h" // For PARAVIEW_USE_MPI

#include "vtkDataArraySelection.h"
#include "vtkDataSetAlgorithm.h"

#ifdef PARAVIEW_USE_MPI
  class vtkMultiProcessController;
#endif

class VTK_EXPORT vtkOGSSpatialStatsFromFile : public vtkRectilinearGridAlgorithm
{
public:
  static vtkOGSSpatialStatsFromFile* New();
  vtkTypeMacro(vtkOGSSpatialStatsFromFile, vtkRectilinearGridAlgorithm);

  // Description:
  // Get the folder where the stats are
  vtkSetStringMacro(FolderName);

  // Description:
  // Get the name of the mask fields to operate
  vtkSetStringMacro(bmask_field);
  vtkSetStringMacro(cmask_field);

  // Description:
  // If true, obtains the statistics per basin and per coast
  vtkGetMacro(per_coast, int);
  vtkSetMacro(per_coast, int);
  vtkBooleanMacro(per_coast, int);

  // Description:
  // The following methods allow selective seleccion of aggregate variables.
  int GetNumberOfStatArrays();
  const char * GetStatArrayName(int index);
  int GetStatArrayIndex(const char* name);
  int GetStatArrayStatus(const char *name);
  void SetStatArrayStatus(const char* name, int status);
  void DisableAllStatArrays();
  void EnableAllStatArrays();

  #ifdef PARAVIEW_USE_MPI
    // Description:
    // Set the controller use in compositing (set to
    // the global controller by default)
    // If not using the default, this must be called before any
    // other methods.
    virtual void SetController(vtkMultiProcessController* controller);
  #endif

protected:
  vtkOGSSpatialStatsFromFile();
  ~vtkOGSSpatialStatsFromFile() override;

  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *) override;

  vtkDataArraySelection* StatDataArraySelection;

  #ifdef PARAVIEW_USE_MPI
    vtkMultiProcessController* Controller;
  #endif

private:
  vtkOGSSpatialStatsFromFile(const vtkOGSSpatialStatsFromFile&) = delete;
  void operator=(const vtkOGSSpatialStatsFromFile&) = delete;

  char *FolderName;
  char *bmask_field, *cmask_field;

  int per_coast, procId, nProcs;
};

#endif
