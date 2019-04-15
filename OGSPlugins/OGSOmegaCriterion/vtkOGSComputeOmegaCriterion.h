// -*- c++ -*-
/*=========================================================================

  Program:   OGSComputeOmegaCriterion
  Module:    vtkOGSComputeOmegaCriterion.h

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef vtkOGSComputeOmegaCriterion_h
#define vtkOGSComputeOmegaCriterion_h

#include "vtkDataSet.h"
#include "vtkRectilinearGridAlgorithm.h"

#include "vtkPVConfig.h" // For PARAVIEW_USE_MPI

#include "../_utils/V3.h"
#include "../_utils/field.h"

#ifdef PARAVIEW_USE_MPI
  class vtkMultiProcessController;
#endif

class VTK_EXPORT vtkOGSComputeOmegaCriterion : public vtkRectilinearGridAlgorithm
{
public:
  static vtkOGSComputeOmegaCriterion* New();
  vtkTypeMacro(vtkOGSComputeOmegaCriterion, vtkRectilinearGridAlgorithm);

  // Description:
  // Get a small number so as not to divide by zero
  vtkSetMacro(epsi, double);
  vtkGetMacro(epsi, double);

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
  vtkOGSComputeOmegaCriterion();
  ~vtkOGSComputeOmegaCriterion() override;

  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *) override;
  int RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;

  #ifdef PARAVIEW_USE_MPI
    vtkMultiProcessController* Controller;
  #endif

private:
  vtkOGSComputeOmegaCriterion(const vtkOGSComputeOmegaCriterion&) = delete;
  void operator=(const vtkOGSComputeOmegaCriterion&) = delete;

  char *field;
  int grad_type, procId, nProcs;
  double epsi;

  bool isReqInfo;

  v3::V3v xyz;               // Stores cell/point coordinates
};

#endif
