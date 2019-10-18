// -*- c++ -*-
/*=========================================================================

  Program:   OGSSelectTools
  Module:    vtkOGSSelectPolygon.h

  Copyright (c) 2019 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef vtkOGSSelectPolygon_h
#define vtkOGSSelectPolygon_h

#include "vtkDataSet.h"
#include "vtkThreshold.h"

#include <string>

#include "vtkPVConfig.h" // For PARAVIEW_USE_MPI

#ifdef PARAVIEW_USE_MPI
  class vtkMultiProcessController;
#endif

class VTK_EXPORT vtkOGSSelectPolygon : public vtkThreshold {
public:
  static vtkOGSSelectPolygon* New();
  vtkTypeMacro(vtkOGSSelectPolygon, vtkThreshold);

  // Description:
  // Get the name of the mask field to operate
  vtkSetStringMacro(strPoints);

  // Description:
  // If false, do not include the meshmask. True by default.
  vtkGetMacro(Invert, int);
  vtkSetMacro(Invert, int);
  vtkBooleanMacro(Invert, int);

  #ifdef PARAVIEW_USE_MPI
    // Description:
    // Set the controller use in compositing (set to
    // the global controller by default)
    // If not using the default, this must be called before any
    // other methods.
    virtual void SetController(vtkMultiProcessController* controller);
  #endif

protected:
  vtkOGSSelectPolygon();
  ~vtkOGSSelectPolygon();

  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *) override;

  #ifdef PARAVIEW_USE_MPI
    vtkMultiProcessController* Controller;
  #endif

private:
  vtkOGSSelectPolygon(const vtkOGSSelectPolygon&) = delete;
  void operator=(const vtkOGSSelectPolygon&) = delete;

  int procId, nProcs, Invert;
  double dfact;
  char *strPoints;
  std::string projName;
};

#endif
