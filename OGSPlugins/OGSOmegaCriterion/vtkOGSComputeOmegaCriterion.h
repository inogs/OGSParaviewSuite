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

#include "../_utils/V3.h"
#include "../_utils/field.h"

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

protected:
  vtkOGSComputeOmegaCriterion();
  ~vtkOGSComputeOmegaCriterion() override;

  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *) override;
  int RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;

private:
  vtkOGSComputeOmegaCriterion(const vtkOGSComputeOmegaCriterion&) = delete;
  void operator=(const vtkOGSComputeOmegaCriterion&) = delete;

  char *field;
  int grad_type;
  double epsi;

  bool isReqInfo;

  v3::V3v xyz;               // Stores cell/point coordinates
};

#endif
