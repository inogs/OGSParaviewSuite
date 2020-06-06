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
  // Get the names of the scalar and vectorial fields to operate
  vtkSetStringMacro(scaFName);
  vtkGetStringMacro(scaFName);
  vtkSetStringMacro(vecFName);
  vtkGetStringMacro(vecFName);

  // Description:
  // If false, do not compute the mask array. False by default.
  vtkGetMacro(computeMask, int);
  vtkSetMacro(computeMask, int);
  vtkBooleanMacro(computeMask, int);

  // Description:
  // Get the name of the mask array
  vtkSetStringMacro(maskName);
  vtkGetStringMacro(maskName);

  // Description:
  // Get the name of the field to computate the mask from
  vtkSetStringMacro(arr2Mask);
  vtkGetStringMacro(arr2Mask);

  // Description:
  // Get the coefficient for the mask
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

  // Description:
  // If true, set as negative the anticyclonic vortices on the mask
  vtkGetMacro(changemask, int);
  vtkSetMacro(changemask, int);
  vtkBooleanMacro(changemask, int);

  // Description:
  // Table settings
  vtkGetMacro(computeBCent, int);
  vtkSetMacro(computeBCent, int);
  vtkBooleanMacro(computeBCent, int);
  vtkGetMacro(computeCent, int);
  vtkSetMacro(computeCent, int);
  vtkBooleanMacro(computeCent, int);
  vtkGetMacro(computeSize, int);
  vtkSetMacro(computeSize, int);
  vtkBooleanMacro(computeSize, int);
  vtkGetMacro(computeRAxis, int);
  vtkSetMacro(computeRAxis, int);
  vtkBooleanMacro(computeRAxis, int);
  vtkGetMacro(computeAStr, int);
  vtkSetMacro(computeAStr, int);
  vtkBooleanMacro(computeAStr, int);
  vtkGetMacro(computeRStr, int);
  vtkSetMacro(computeRStr, int);
  vtkBooleanMacro(computeRStr, int);
  vtkGetMacro(computeCirc, int);
  vtkSetMacro(computeCirc, int);
  vtkBooleanMacro(computeCirc, int);

  // Description:
  // Get the names of the normal array
//  vtkSetStringMacro(normName);
//  vtkGetStringMacro(normName);

  // Let the user select a multiplier factor for the depth
  vtkGetMacro(dfact, double);
  vtkSetMacro(dfact, double);

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
  int computeBCent, computeCent, computeSize, computeRAxis, computeAStr, computeRStr, computeCirc;
  int procId, nProcs;

  char *scaFName, *vecFName, *arr2Mask, *maskName; //*normName;

  std::string projName;

  double coef, dfact;
};

#endif