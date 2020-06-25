/*=========================================================================

  Program:   OGSClimatology
  Module:    vtkOGSClimatology.h

  Copyright (c) 2020 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef vtkOGSClimatology_h
#define vtkOGSClimatology_h

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

class VTK_EXPORT vtkOGSClimatology : public vtkDataSetAlgorithm
{
public:
  static vtkOGSClimatology* New();
  vtkTypeMacro(vtkOGSClimatology, vtkDataSetAlgorithm);

  // Description:
  // Selection of the requestor
  vtkGetMacro(ReqType, int);
  vtkSetMacro(ReqType, int);

  // Description:
  // Selection of the algorithm
  vtkGetMacro(use_files, int);
  vtkSetMacro(use_files, int);
  vtkBooleanMacro(use_files, int);

  vtkStringArray *GetTimeValues();

  #ifdef PARAVIEW_USE_MPI
  vtkGetObjectMacro(Controller, vtkMultiProcessController);
  virtual void SetController(vtkMultiProcessController*);
  #endif

protected:
  vtkOGSClimatology();
  ~vtkOGSClimatology();

  int RequestInformation(vtkInformation *, vtkInformationVector **, vtkInformationVector *) override;
  int RequestUpdateExtent(vtkInformation *, vtkInformationVector **, vtkInformationVector *) override;
  int RequestData(vtkInformation *, vtkInformationVector **,vtkInformationVector *) override;

  #ifdef PARAVIEW_USE_MPI
  vtkMultiProcessController* Controller;
  #endif

private:
  vtkOGSClimatology(const vtkOGSClimatology&) = delete;
  void operator=(const vtkOGSClimatology&) = delete;

  int PipelineIterationAlgorithm(vtkInformation *, vtkDataSet *, vtkDataSet *);
  int FileIterationAlgorithm(vtkInformation *, vtkDataSet *, vtkDataSet *);

  Time::TimeList TL;      // TimeList containing all the instants

  std::vector<int> instants;   // Instant ID to loop
  std::vector<double> weights; // Weights for the instants

  double sum_weights;
  bool TL_computed, use_files;

  int procId, nProcs, ReqType, CurrentTimeIndex;
};

#endif
