/*=========================================================================

  Program:   OGSVariableAggregator
  Module:    vtkOGSVariableAggregator.h

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef vtkOGSVariableAggregator_h
#define vtkOGSVariableAggregator_h

#include "vtkDataArraySelection.h"
#include "vtkDataSetAlgorithm.h"

class VTK_EXPORT vtkOGSVariableAggregator : public vtkDataSetAlgorithm
{
public:
  static vtkOGSVariableAggregator* New();
  vtkTypeMacro(vtkOGSVariableAggregator, vtkDataSetAlgorithm);

  // Description:
  // If false, the aggregated variables will not be deleted
  vtkGetMacro(deleteVars, int);
  vtkSetMacro(deleteVars, int);
  vtkBooleanMacro(deleteVars, int);

  // Description:
  // The following methods allow selective seleccion of aggregate variables.
  int GetNumberOfVarArrays();
  const char * GetVarArrayName(int index);
  int GetVarArrayIndex(const char* name);
  int GetVarArrayStatus(const char *name);
  void SetVarArrayStatus(const char* name, int status);
  void DisableAllVarArrays();
  void EnableAllVarArrays();

protected:
  vtkOGSVariableAggregator();
  ~vtkOGSVariableAggregator();

  int RequestData(vtkInformation *, vtkInformationVector **,vtkInformationVector *) override;

  vtkDataArraySelection* VarDataArraySelection; // Stores name of the aggregated variable
  vtkDataArraySelection* AgrDataArraySelection; // Stores the variables to aggregate separated by ";"

private:
  vtkOGSVariableAggregator(const vtkOGSVariableAggregator&) = delete;
  void operator=(const vtkOGSVariableAggregator&) = delete;

  class vtkVectorOfArrays;

  int deleteVars;
};

#endif
