/*=========================================================================

  Program:   OGSTimeAggregator
  Module:    vtkOGSTimeAggregator.h

  Copyright (c) 2020 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef vtkOGSTimeAggregator_h
#define vtkOGSTimeAggregator_h

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

class VTK_EXPORT vtkOGSTimeAggregator : public vtkDataSetAlgorithm
{
public:
  static vtkOGSTimeAggregator* New();
  vtkTypeMacro(vtkOGSTimeAggregator, vtkDataSetAlgorithm);

  // Description:
  // Selection of the requestor
  vtkGetMacro(ReqType, int);
  vtkSetMacro(ReqType, int);

  // Description:
  // Selection of the weekday for the weekly list
  vtkGetMacro(weekday, int);
  vtkSetMacro(weekday, int);

  // Description:
  // Selection of the algorithm
  vtkGetMacro(use_files, int);
  vtkSetMacro(use_files, int);
  vtkBooleanMacro(use_files, int);

  #ifdef PARAVIEW_USE_MPI
  vtkGetObjectMacro(Controller, vtkMultiProcessController);
  virtual void SetController(vtkMultiProcessController*);
  #endif

protected:
  vtkOGSTimeAggregator();
  ~vtkOGSTimeAggregator();

  int RequestInformation(vtkInformation *, vtkInformationVector **, vtkInformationVector *) override;
  int RequestUpdateExtent(vtkInformation *, vtkInformationVector **, vtkInformationVector *) override;
  int RequestData(vtkInformation *, vtkInformationVector **,vtkInformationVector *) override;

  #ifdef PARAVIEW_USE_MPI
  vtkMultiProcessController* Controller;
  #endif

private:
  vtkOGSTimeAggregator(const vtkOGSTimeAggregator&) = delete;
  void operator=(const vtkOGSTimeAggregator&) = delete;

  int PipelineIterationAlgorithm(vtkInformation *, vtkDataSet *, vtkDataSet *);
  int FileIterationAlgorithm(vtkInformation *, vtkDataSet *, vtkDataSet *);

  Time::TimeList TL;      // TimeList containing all the instants
  Time::REQ_LIST ReqList; // Requestor list


  std::vector<int> instants;   // Instant ID to loop
  std::vector<double> weights; // Weights for the instants

  double sum_weights;
  bool TL_computed, use_files;

  int procId, nProcs, weekday, ReqType, CurrentTimeIndex;
};

#endif
