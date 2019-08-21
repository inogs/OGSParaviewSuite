// -*- c++ -*-
/*=========================================================================

  Program:   OGSSelectTools
  Module:    vtkOGSSelectCoast.h

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef vtkOGSSelectCoast_h
#define vtkOGSSelectCoast_h

#include "vtkDataSet.h"
#include "vtkThreshold.h"

#include "vtkPVConfig.h" // For PARAVIEW_USE_MPI

class vtkDataArraySelection;

#ifdef PARAVIEW_USE_MPI
  class vtkMultiProcessController;
#endif

class VTK_EXPORT vtkOGSSelectCoast : public vtkThreshold {
public:
  static vtkOGSSelectCoast* New();
  vtkTypeMacro(vtkOGSSelectCoast, vtkThreshold);

  void PrintSelf(ostream& os, vtkIndent indent) VTK_OVERRIDE;

  // Description:
  // Get the name of the mask field to operate
  vtkSetStringMacro(mask_field);

  // Description:
  // The following methods allow selective seleccion of the basins.
  int GetNumberOfCoastsArrays();
  const char * GetCoastsArrayName(int index);
  int GetCoastsArrayIndex(const char* name);
  int GetCoastsArrayStatus(const char *name);
  void SetCoastsArrayStatus(const char* name, int status);
  void DisableAllCoastsArrays();
  void EnableAllCoastsArrays();

  #ifdef PARAVIEW_USE_MPI
    // Description:
    // Set the controller use in compositing (set to
    // the global controller by default)
    // If not using the default, this must be called before any
    // other methods.
    virtual void SetController(vtkMultiProcessController* controller);
  #endif

protected:
  vtkOGSSelectCoast();
  ~vtkOGSSelectCoast();

  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *) override;

  vtkDataArraySelection* CoastsDataArraySelection;

  #ifdef PARAVIEW_USE_MPI
    vtkMultiProcessController* Controller;
  #endif

private:
  vtkOGSSelectCoast(const vtkOGSSelectCoast&) = delete;
  void operator=(const vtkOGSSelectCoast&) = delete;

  int procId, nProcs;
  char *mask_field;
};

#endif
