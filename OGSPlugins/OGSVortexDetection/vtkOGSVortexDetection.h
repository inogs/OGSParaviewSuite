// -*- c++ -*-
/*=========================================================================

  Program:   OGSVortexDetection
  Module:    vtkOGSVortexDetection.h

  Copyright (c) 2020 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef vtkOGSVortexDetection_h
#define vtkOGSVortexDetection_h

#include "vtkDataSet.h"
#include "vtkFiltersGeneralModule.h" // For export macro
#include "vtkDataSetAlgorithm.h"
#include "vtkTableAlgorithm.h"

#include "vtkPVConfig.h" // For PARAVIEW_USE_MPI

#include <string>

#include "vortex.h"

class vtkStringArray;

#ifdef PARAVIEW_USE_MPI
  class vtkMultiProcessController;
#endif

class vtkOGSVortexDetection : public vtkDataSetAlgorithm
{
public:
  static vtkOGSVortexDetection* New();
  vtkTypeMacro(vtkOGSVortexDetection, vtkDataSetAlgorithm);
  
  // Description:
  // Get the name of the mask field to operate
  vtkSetStringMacro(arrName);

  // Description:
  // Get the name of the mask field to operate
  vtkSetStringMacro(omegArrName);

  // Description:
  // Get the name of the mask field to operate
  vtkSetStringMacro(maskArrName1);
  vtkSetStringMacro(maskArrName2);

  // Description:
  // If false, do not compute the mask array. False by default.
  vtkGetMacro(computeMask, int);
  vtkSetMacro(computeMask, int);
  vtkBooleanMacro(computeMask, int);

  // Description:
  // Get the coefficient for the Okubo-Weiss mask
  vtkSetMacro(coef, double);
  vtkGetMacro(coef, double);

  // Description:
  // Control the maximum number of repetitions for 
  // the vortex detection algorithm
  vtkSetMacro(maxreps, int);
  vtkGetMacro(maxreps, int);

  // Description:
  // Control the maximum number of iterations for 
  // the vortex detection algorithm
  vtkSetMacro(maxiter, int);
  vtkGetMacro(maxiter, int);

  // Description:
  // Minimum number of elements to be considered a vortex
  vtkSetMacro(minres, int);
  vtkGetMacro(minres, int);

  // Let the user select a multiplier factor for the depth
  vtkGetMacro(dfact, double);
  vtkSetMacro(dfact, double);

  // Description:
  // If true, set as negative the anticyclonic vortices on the mask
  vtkGetMacro(changemask, int);
  vtkSetMacro(changemask, int);
  vtkBooleanMacro(changemask, int);

  // Description:
  // Get the output data object for a port on this algorithm.
  vtkDataSet* GetOutputA();
  vtkTable*   GetOutputB();
  
  void SetOutput(vtkDataObject* d);

  vtkStringArray * GetProjections();
  void SetProjection(const char *proj);

  #ifdef PARAVIEW_USE_MPI
    // Description:
    // Set the controller use in compositing (set to
    // the global controller by default)
    // If not using the default, this must be called before any
    // other methods.
    virtual void SetController(vtkMultiProcessController* controller);
  #endif

protected:
  vtkOGSVortexDetection();
  ~vtkOGSVortexDetection() override;

  int RequestDataObject(vtkInformation*,vtkInformationVector**,vtkInformationVector*);
  int RequestUpdateExtent(vtkInformation*,vtkInformationVector**,vtkInformationVector*);
  int FillOutputPortInformation(int, vtkInformation*);

  int RequestInformation(vtkInformation*,vtkInformationVector**,vtkInformationVector*);
  int RequestData(vtkInformation*,vtkInformationVector**,vtkInformationVector*);

  vtkStringArray *Projections;

  #ifdef PARAVIEW_USE_MPI
    vtkMultiProcessController* Controller;
  #endif

private:
  vtkOGSVortexDetection(const vtkOGSVortexDetection &) = delete;
  void operator=(const vtkOGSVortexDetection &) = delete;

  int computeMaskArray(bool,vtkDataSet*);
  int computeDetectionMask(bool,vtkDataSet*,vortex::VortexList&);
  int computeVortexProperties(bool,vortex::VortexList&,vtkDataSet*,vtkTable*);

  int maxreps, maxiter, minres, computeMask, changemask;
  int procId, nProcs;

  char *arrName, *omegArrName, *maskArrName1, *maskArrName2;

  std::string projName;

  double coef, dfact;
};

#endif