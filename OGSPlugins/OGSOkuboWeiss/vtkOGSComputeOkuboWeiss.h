// -*- c++ -*-
/*=========================================================================

  Program:   OGSComputeOkuboWeiss
  Module:    vtkOGSComputeOkuboWeiss.h

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef vtkOGSComputeOkuboWeiss_h
#define vtkOGSComputeOkuboWeiss_h

#include "vtkDataSet.h"
#include "vtkRectilinearGridAlgorithm.h"

#include "vtkPVConfig.h" // For PARAVIEW_USE_MPI

#include "../_utils/V3.h"
#include "../_utils/field.h"

#ifdef PARAVIEW_USE_MPI
  class vtkMultiProcessController;
#endif

class VTK_EXPORT vtkOGSComputeOkuboWeiss : public vtkRectilinearGridAlgorithm
{
public:
  static vtkOGSComputeOkuboWeiss* New();
  vtkTypeMacro(vtkOGSComputeOkuboWeiss, vtkRectilinearGridAlgorithm);

  // Description:
  // Get the coefficient for the Okubo-Weiss mask
  vtkSetMacro(coef, double);
  vtkGetMacro(coef, double);

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
  vtkOGSComputeOkuboWeiss();
  ~vtkOGSComputeOkuboWeiss() override;

  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *) override;
  int RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;

  #ifdef PARAVIEW_USE_MPI
    vtkMultiProcessController* Controller;
  #endif

private:
  vtkOGSComputeOkuboWeiss(const vtkOGSComputeOkuboWeiss&) = delete;
  void operator=(const vtkOGSComputeOkuboWeiss&) = delete;

  char *field;
  double coef;
  int grad_type, procId, nProcs;

  bool isReqInfo;

  v3::V3v xyz;               // Stores cell/point coordinates
};

#endif
