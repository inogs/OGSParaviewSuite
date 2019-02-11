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

protected:
  vtkOGSComputeLambda2Criterion();
  ~vtkOGSComputeLambda2Criterion() override;

  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *) override;

private:
  vtkOGSComputeLambda2Criterion(const vtkOGSComputeLambda2Criterion&) = delete;
  void operator=(const vtkOGSComputeLambda2Criterion&) = delete;

  char *field;
  int grad_type;
};

#endif
