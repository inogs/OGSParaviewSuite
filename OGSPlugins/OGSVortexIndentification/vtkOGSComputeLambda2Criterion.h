// -*- c++ -*-
/*=========================================================================

  Program:   OGSComputeLambda2Criterion
  Module:    vtkOGSComputeLambda2Criterion.h

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef vtkOGSComputeQCriterion_h
#define vtkOGSComputeQCriterion_h

#include "vtkDataSet.h"
#include "vtkRectilinearGridAlgorithm.h"

#include "vtkPVConfig.h" // For PARAVIEW_USE_MPI

#include "V3.h"
#include "field.h"

#ifdef PARAVIEW_USE_MPI
  class vtkMultiProcessController;
#endif

class VTK_EXPORT vtkOGSComputeLambda2Criterion : public vtkRectilinearGridAlgorithm
{
public:
  static vtkOGSComputeLambda2Criterion* New();
  vtkTypeMacro(vtkOGSComputeLambda2Criterion, vtkRectilinearGridAlgorithm);

  // Description:
  // Gradient method for computing Okubo-Weiss
  vtkSetMacro(grad_type, int);
  vtkGetMacro(grad_type, int);

  // Description:
  // Get the name of the velocity field
  vtkSetStringMacro(field);

  #ifdef PARAVIEW_USE_MPI
    // Description:
    // Set the controller use in compositing (set to
    // the global controller by default)
    // If not using the default, this must be called before any
    // other methods.
    virtual void SetController(vtkMultiProcessController* controller);
  #endif

protected:
  vtkOGSComputeLambda2Criterion();
  ~vtkOGSComputeLambda2Criterion() override;

  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *) override;
  int RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;

  #ifdef PARAVIEW_USE_MPI
    vtkMultiProcessController* Controller;
  #endif

private:
  vtkOGSComputeLambda2Criterion(const vtkOGSComputeLambda2Criterion&) = delete;
  void operator=(const vtkOGSComputeLambda2Criterion&) = delete;

  char *field;
  int grad_type, procId, nProcs;

  bool isReqInfo;

  v3::V3v xyz;               // Stores cell/point coordinates
};

#endif
