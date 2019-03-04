// -*- c++ -*-
/*=========================================================================

  Program:   OGSComputeQCriterion
  Module:    vtkOGSComputeQCriterion.h

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

#include "../_utils/V3.h"
#include "../_utils/field.h"

class VTK_EXPORT vtkOGSComputeQCriterion : public vtkRectilinearGridAlgorithm
{
public:
  static vtkOGSComputeQCriterion* New();
  vtkTypeMacro(vtkOGSComputeQCriterion, vtkRectilinearGridAlgorithm);

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

protected:
  vtkOGSComputeQCriterion();
  ~vtkOGSComputeQCriterion() override;

  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *) override;
  int RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;

private:
  vtkOGSComputeQCriterion(const vtkOGSComputeQCriterion&) = delete;
  void operator=(const vtkOGSComputeQCriterion&) = delete;

  char *field;
  double coef;
  int grad_type;

  bool isReqInfo;

  v3::V3v xyz;               // Stores cell/point coordinates
};

#endif
