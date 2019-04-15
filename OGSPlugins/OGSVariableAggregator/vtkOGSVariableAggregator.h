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

#include "vtkPVConfig.h" // For PARAVIEW_USE_MPI

#include <map>

#ifdef PARAVIEW_USE_MPI
  class vtkMultiProcessController;
#endif

class VTK_EXPORT vtkOGSVariableAggregator : public vtkDataSetAlgorithm
{
public:
  static vtkOGSVariableAggregator* New();
  vtkTypeMacro(vtkOGSVariableAggregator, vtkDataSetAlgorithm);

  // Description:
  // Get the name of the XML file to read
  vtkSetStringMacro(FileName);

  // Description:
  // If false, the aggregated variables will not be deleted
  vtkGetMacro(deleteVars, int);
  vtkSetMacro(deleteVars, int);
  vtkBooleanMacro(deleteVars, int);

  // Description:
  // Parse XML text in TextBox
  vtkSetStringMacro(XMLText);

  // Description:
  // The following methods allow selective seleccion of aggregate variables.
  int GetNumberOfVarArrays();
  const char * GetVarArrayName(int index);
  int GetVarArrayIndex(const char* name);
  int GetVarArrayStatus(const char *name);
  void SetVarArrayStatus(const char* name, int status);
  void DisableAllVarArrays();
  void EnableAllVarArrays();

  #ifdef PARAVIEW_USE_MPI
    // Description:
    // Set the controller use in compositing (set to
    // the global controller by default)
    // If not using the default, this must be called before any
    // other methods.
    virtual void SetController(vtkMultiProcessController* controller);
  #endif

protected:
  vtkOGSVariableAggregator();
  ~vtkOGSVariableAggregator();

  int RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int RequestData(vtkInformation *, vtkInformationVector **,vtkInformationVector *) override;

  vtkDataArraySelection* VarDataArraySelection; // Stores name of the aggregated variable
  std::map<std::string, std::string> AggrVar;   // Stores the variables to aggregate separated by ";"

  #ifdef PARAVIEW_USE_MPI
    vtkMultiProcessController* Controller;
  #endif

private:
  vtkOGSVariableAggregator(const vtkOGSVariableAggregator&) = delete;
  void operator=(const vtkOGSVariableAggregator&) = delete;

  void ParseXML();
  void SetAggrVarsText();

  int deleteVars, procId, nProcs;
  char *FileName, *XMLText;
};

#endif
