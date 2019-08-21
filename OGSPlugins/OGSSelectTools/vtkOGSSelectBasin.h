// -*- c++ -*-
/*=========================================================================

  Program:   OGSSelectTools
  Module:    vtkOGSSelectBasin.h

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef vtkOGSSelectBasin_h
#define vtkOGSSelectBasin_h

#include "vtkDataSet.h"
#include "vtkThreshold.h"

#include "vtkPVConfig.h" // For PARAVIEW_USE_MPI

class vtkDataArraySelection;

#ifdef PARAVIEW_USE_MPI
  class vtkMultiProcessController;
#endif

class VTK_EXPORT vtkOGSSelectBasin : public vtkThreshold {
public:
  static vtkOGSSelectBasin* New();
  vtkTypeMacro(vtkOGSSelectBasin, vtkThreshold);

  void PrintSelf(ostream& os, vtkIndent indent) VTK_OVERRIDE;

  // Description:
  // Get the name of the mask field to operate
  vtkSetStringMacro(mask_field);

  // Description:
  // The following methods allow selective seleccion of the basins.
  int GetNumberOfBasinsArrays();
  const char * GetBasinsArrayName(int index);
  int GetBasinsArrayIndex(const char* name);
  int GetBasinsArrayStatus(const char *name);
  void SetBasinsArrayStatus(const char* name, int status);
  void DisableAllBasinsArrays();
  void EnableAllBasinsArrays();

  #ifdef PARAVIEW_USE_MPI
    // Description:
    // Set the controller use in compositing (set to
    // the global controller by default)
    // If not using the default, this must be called before any
    // other methods.
    virtual void SetController(vtkMultiProcessController* controller);
  #endif

protected:
  vtkOGSSelectBasin();
  ~vtkOGSSelectBasin() override;

  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *) override;

  vtkDataArraySelection* BasinsDataArraySelection;

  #ifdef PARAVIEW_USE_MPI
    vtkMultiProcessController* Controller;
  #endif

private:
  vtkOGSSelectBasin(const vtkOGSSelectBasin&) = delete;
  void operator=(const vtkOGSSelectBasin&) = delete;

  int procId, nProcs;
  char *mask_field;
};

#endif
