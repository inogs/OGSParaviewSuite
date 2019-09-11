/*=========================================================================

  Program:   OGSTimeStatistics
  Module:    vtkOGSTimeStatistics.h

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef vtkOGSTimeStatistics_h
#define vtkOGSTimeStatistics_h

#include "vtkDataSetAlgorithm.h"

#include <string>

#include "vtkPVConfig.h" // For PARAVIEW_USE_MPI

class vtkStringArray;

#ifdef PARAVIEW_USE_MPI
class vtkMultiProcessController;
#endif

class VTK_EXPORT vtkOGSTimeStatistics : public vtkDataSetAlgorithm
{
public:
  static vtkOGSTimeStatistics* New();
  vtkTypeMacro(vtkOGSTimeStatistics, vtkDataSetAlgorithm);

  void SetStartTime(const char *tstep);
  void SetEndTime(const char *tstep);

  vtkStringArray *GetTimeValues();

  #ifdef PARAVIEW_USE_MPI
  vtkGetObjectMacro(Controller, vtkMultiProcessController);
  virtual void SetController(vtkMultiProcessController*);
  #endif

protected:
  vtkOGSTimeStatistics();
  ~vtkOGSTimeStatistics();

  int RequestInformation(vtkInformation *, vtkInformationVector **, vtkInformationVector *) override;
  int RequestData(vtkInformation *, vtkInformationVector **,vtkInformationVector *) override;

  vtkStringArray *TimeValues;

  #ifdef PARAVIEW_USE_MPI
  vtkMultiProcessController* Controller;
  #endif

private:
  vtkOGSTimeStatistics(const vtkOGSTimeStatistics&) = delete;
  void operator=(const vtkOGSTimeStatistics&) = delete;

  int procId, nProcs;
  int ii_start, ii_end;

  std::string tstep_st, tstep_ed;
};

#endif
