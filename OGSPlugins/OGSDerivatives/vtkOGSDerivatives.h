// -*- c++ -*-
/*=========================================================================

  Program:   OGSDerivatives
  Module:    vtkOGSDerivatives.h

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef vtkOGSDerivatives_h
#define vtkOGSDerivatives_h

#include "vtkDataSet.h"
#include "vtkFloatArray.h"
#include "vtkRectilinearGridAlgorithm.h"

#include "../_utils/V3.h"
#include "../_utils/field.h"

class VTK_EXPORT vtkOGSDerivatives : public vtkRectilinearGridAlgorithm
{
public:
  static vtkOGSDerivatives* New();
  vtkTypeMacro(vtkOGSDerivatives, vtkRectilinearGridAlgorithm);

  // Description:
  // Gradient method for computing Okubo-Weiss
  vtkSetMacro(grad_type, int);
  vtkGetMacro(grad_type, int);

  // Description:
  // Computation of the divergence
  vtkGetMacro(ComputeDivergence, int);
  vtkSetMacro(ComputeDivergence, int);
  vtkBooleanMacro(ComputeDivergence, int);

  // Description:
  // Computation of the curl
  vtkGetMacro(ComputeCurl, int);
  vtkSetMacro(ComputeCurl, int);
  vtkBooleanMacro(ComputeCurl, int);

  // Description:
  // Computation of the Q-criterion
  vtkGetMacro(ComputeQ, int);
  vtkSetMacro(ComputeQ, int);
  vtkBooleanMacro(ComputeQ, int);

  // Description:
  // Get the name of the velocity field
  vtkSetStringMacro(field);

protected:
  vtkOGSDerivatives();
  ~vtkOGSDerivatives() override;

  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *) override;
  int RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;

private:
  vtkOGSDerivatives(const vtkOGSDerivatives&) = delete;
  void operator=(const vtkOGSDerivatives&) = delete;

  char *field;

  int grad_type;
  
  bool isReqInfo;
  bool ComputeDivergence, ComputeCurl, ComputeQ;

  v3::V3v xyz;               // Stores cell/point coordinates
};

#endif
