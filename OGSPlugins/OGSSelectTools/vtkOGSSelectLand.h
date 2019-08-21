// -*- c++ -*-
/*=========================================================================

  Program:   OGSSelectTools
  Module:    vtkOGSSelectLand.h

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef vtkOGSSelectLand_h
#define vtkOGSSelectLand_h

#include "vtkDataSet.h"
#include "vtkThreshold.h"

#include "vtkPVConfig.h" // For PARAVIEW_USE_MPI

#ifdef PARAVIEW_USE_MPI
  class vtkMultiProcessController;
#endif

class VTK_EXPORT vtkOGSSelectLand : public vtkThreshold {
public:
  static vtkOGSSelectLand* New();
  vtkTypeMacro(vtkOGSSelectLand, vtkThreshold);

  void PrintSelf(ostream& os, vtkIndent indent) VTK_OVERRIDE;

  // Description:
  // Get the name of the mask field to operate
  vtkSetStringMacro(mask_field);

  #ifdef PARAVIEW_USE_MPI
    // Description:
    // Set the controller use in compositing (set to
    // the global controller by default)
    // If not using the default, this must be called before any
    // other methods.
    virtual void SetController(vtkMultiProcessController* controller);
  #endif

protected:
  vtkOGSSelectLand();
  ~vtkOGSSelectLand();

  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *) override;

  #ifdef PARAVIEW_USE_MPI
    vtkMultiProcessController* Controller;
  #endif

private:
  vtkOGSSelectLand(const vtkOGSSelectLand&) = delete;
  void operator=(const vtkOGSSelectLand&) = delete;

  int procId, nProcs;
  char *mask_field;
};

#endif
