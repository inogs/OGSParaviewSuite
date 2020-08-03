/*=========================================================================

  Program:   OGSSelectTimePeriod
  Module:    vtkOGSSelectTimePeriod.h

  Copyright (c) 2020 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef vtkOGSSelectTimePeriod_h
#define vtkOGSSelectTimePeriod_h

#include "vtkDataSetAlgorithm.h"
#include "vtkStringArray.h"

#include <string>
#include <vector>
#include "TimeInterval.h"
#include "TimeList.h"

#include "vtkPVConfig.h" // For PARAVIEW_USE_MPI

class vtkStringArray;

#ifdef PARAVIEW_USE_MPI
class vtkMultiProcessController;
#endif

class VTK_EXPORT vtkOGSSelectTimePeriod : public vtkDataSetAlgorithm
{
public:
  static vtkOGSSelectTimePeriod* New();
  vtkTypeMacro(vtkOGSSelectTimePeriod, vtkDataSetAlgorithm);

  void SetStartTI(const char *tstep);
  void SetEndTI(const char *tstep);

  vtkStringArray *GetTimeValues();

  #ifdef PARAVIEW_USE_MPI
  vtkGetObjectMacro(Controller, vtkMultiProcessController);
  virtual void SetController(vtkMultiProcessController*);
  #endif

protected:
  vtkOGSSelectTimePeriod();
  ~vtkOGSSelectTimePeriod();

  int RequestInformation(vtkInformation *, vtkInformationVector **, vtkInformationVector *) override;
  int RequestData(vtkInformation *, vtkInformationVector **,vtkInformationVector *) override;

  #ifdef PARAVIEW_USE_MPI
  vtkMultiProcessController* Controller;
  #endif

private:
  vtkOGSSelectTimePeriod(const vtkOGSSelectTimePeriod&) = delete;
  void operator=(const vtkOGSSelectTimePeriod&) = delete;

  Time::TimeInterval TI; // TimeInterval for the generic requestor
  Time::TimeList TL;     // TimeList containing all the instants

  std::vector<int> instants;   // Instant ID to loop
  std::vector<double> weights; // Weights for the instants

  bool TL_computed;

  int procId, nProcs, CurrentTimeIndex;
};

#endif
