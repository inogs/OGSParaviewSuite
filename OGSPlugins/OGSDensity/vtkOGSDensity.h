/*=========================================================================

  Program:   OGSDensity
  Module:    vtkOGSDensity.h

  Copyright (c) 2020 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef vtkOGSDensity_h
#define vtkOGSDensity_h

#include "vtkDataSetAlgorithm.h"

#include "vtkPVConfig.h" // For PARAVIEW_USE_MPI

#ifdef PARAVIEW_USE_MPI
  class vtkMultiProcessController;
#endif

class VTK_EXPORT vtkOGSDensity : public vtkDataSetAlgorithm
{
public:
  static vtkOGSDensity* New();
  vtkTypeMacro(vtkOGSDensity, vtkDataSetAlgorithm);

  // Description:
  // Method of computation
  vtkGetMacro(method, int);
  vtkSetMacro(method, int); 

  // Description:
  // Parse temperature and salinity array names
  vtkSetStringMacro(Tarrname);
  vtkSetStringMacro(Sarrname);

  #ifdef PARAVIEW_USE_MPI
    // Description:
    // Set the controller use in compositing (set to
    // the global controller by default)
    // If not using the default, this must be called before any
    // other methods.
    virtual void SetController(vtkMultiProcessController* controller);
  #endif

protected:
  vtkOGSDensity();
  ~vtkOGSDensity();

  int RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int RequestData(vtkInformation *, vtkInformationVector **,vtkInformationVector *) override;

  #ifdef PARAVIEW_USE_MPI
    vtkMultiProcessController* Controller;
  #endif

private:
  vtkOGSDensity(const vtkOGSDensity&) = delete;
  void operator=(const vtkOGSDensity&) = delete;

  int method, procId, nProcs;
  char *Tarrname, *Sarrname;
};

#endif
