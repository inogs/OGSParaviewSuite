/*=========================================================================

  Program:   OGSAnnotateDateTime
  Module:    vtkOGSAnnotateDateTime.h

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef vtkOGSAnnotateDateTime_h
#define vtkOGSAnnotateDateTime_h

#include "vtkTimeToTextConvertor.h"

#ifdef PARAVIEW_USE_MPI
  class vtkMultiProcessController;
#endif

class VTK_EXPORT vtkOGSAnnotateDateTime : public vtkTimeToTextConvertor
{
public:
  static vtkOGSAnnotateDateTime* New();
  vtkTypeMacro(vtkOGSAnnotateDateTime, vtkTimeToTextConvertor);

  vtkSetStringMacro(TimeFormat);
  vtkGetStringMacro(TimeFormat);

  vtkGetMacro(useMetadata, int);
  vtkSetMacro(useMetadata, int);
  vtkBooleanMacro(useMetadata, int);

  #ifdef PARAVIEW_USE_MPI
    // Description:
    // Set the controller use in compositing (set to
    // the global controller by default)
    // If not using the default, this must be called before any
    // other methods.
    virtual void SetController(vtkMultiProcessController* controller);
  #endif

protected:
  vtkOGSAnnotateDateTime();
  ~vtkOGSAnnotateDateTime();

  int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;

  #ifdef PARAVIEW_USE_MPI
    vtkMultiProcessController* Controller;
  #endif

private:
  vtkOGSAnnotateDateTime(const vtkOGSAnnotateDateTime&) = delete;
  void operator=(const vtkOGSAnnotateDateTime&) = delete;

  char *TimeFormat;
  int useMetadata, procId, nProcs;
};

#endif
