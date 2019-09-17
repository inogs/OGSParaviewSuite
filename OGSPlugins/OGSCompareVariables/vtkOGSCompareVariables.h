/*=========================================================================

  Program:   OGSCompareVariables
  Module:    vtkOGSCompareVariables.h

  Copyright (c) 2019 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef vtkOGSCompareVariables_h
#define vtkOGSCompareVariables_h

#include "vtkDataArraySelection.h"
#include "vtkDataSetAlgorithm.h"

#include "vtkPVConfig.h" // For PARAVIEW_USE_MPI

#ifdef PARAVIEW_USE_MPI
  class vtkMultiProcessController;
#endif

class VTK_EXPORT vtkOGSCompareVariables : public vtkDataSetAlgorithm
{
public:
  static vtkOGSCompareVariables* New();
  vtkTypeMacro(vtkOGSCompareVariables, vtkDataSetAlgorithm);

  // Description:
  // Mode of comparison
  vtkGetMacro(mode, int);
  vtkSetMacro(mode, int);

  // Description:
  // Mode of comparison
  vtkGetMacro(defval, double);
  vtkSetMacro(defval, double);  

  // Description:
  // Parse variable text in TextBox
  vtkSetStringMacro(variables);

  #ifdef PARAVIEW_USE_MPI
    // Description:
    // Set the controller use in compositing (set to
    // the global controller by default)
    // If not using the default, this must be called before any
    // other methods.
    virtual void SetController(vtkMultiProcessController* controller);
  #endif

protected:
  vtkOGSCompareVariables();
  ~vtkOGSCompareVariables();

  int RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int RequestData(vtkInformation *, vtkInformationVector **,vtkInformationVector *) override;

  #ifdef PARAVIEW_USE_MPI
    vtkMultiProcessController* Controller;
  #endif

private:
  vtkOGSCompareVariables(const vtkOGSCompareVariables&) = delete;
  void operator=(const vtkOGSCompareVariables&) = delete;

  int mode, procId, nProcs;
  char *variables;
  double defval;
};

#endif
