// -*- c++ -*-
/*=========================================================================

  Program:   OGSSelectOkuboWeiss
  Module:    vtkOGSSelectOkuboWeiss.h

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef vtkOGSSelectOkuboWeiss_h
#define vtkOGSSelectOkuboWeiss_h

#include "vtkDataSet.h"
#include "vtkThreshold.h"

#include "vtkPVConfig.h" // For PARAVIEW_USE_MPI

class vtkDataArraySelection;

#ifdef PARAVIEW_USE_MPI
  class vtkMultiProcessController;
#endif

class VTK_EXPORT vtkOGSSelectOkuboWeiss : public vtkThreshold {
public:
  static vtkOGSSelectOkuboWeiss* New();
  vtkTypeMacro(vtkOGSSelectOkuboWeiss, vtkThreshold);

  void PrintSelf(ostream& os, vtkIndent indent) VTK_OVERRIDE;

  // Description:
  // Get the name of the mask field to operate
  vtkSetStringMacro(mask_field);

  // Description:
  // The following methods allow selective seleccion of the basins.
  int GetNumberOfOWArrays();
  const char * GetOWArrayName(int index);
  int GetOWArrayIndex(const char* name);
  int GetOWArrayStatus(const char *name);
  void SetOWArrayStatus(const char* name, int status);
  void DisableAllOWArrays();
  void EnableAllOWArrays();

  #ifdef PARAVIEW_USE_MPI
    // Description:
    // Set the controller use in compositing (set to
    // the global controller by default)
    // If not using the default, this must be called before any
    // other methods.
    virtual void SetController(vtkMultiProcessController* controller);
  #endif

protected:
  vtkOGSSelectOkuboWeiss();
  ~vtkOGSSelectOkuboWeiss() override;

  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *) override;

  vtkDataArraySelection* OWDataArraySelection;

  #ifdef PARAVIEW_USE_MPI
    vtkMultiProcessController* Controller;
  #endif

private:
  vtkOGSSelectOkuboWeiss(const vtkOGSSelectOkuboWeiss&) = delete;
  void operator=(const vtkOGSSelectOkuboWeiss&) = delete;

  int procId, nProcs;
  char *mask_field;
};

#endif
