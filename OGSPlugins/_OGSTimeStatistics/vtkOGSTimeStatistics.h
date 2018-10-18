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

#include "vtkStringArray.h"
#include "vtkDataSetAlgorithm.h"

class VTK_EXPORT vtkOGSTimeStatistics : public vtkDataSetAlgorithm
{
public:
  static vtkOGSTimeStatistics* New();
  vtkTypeMacro(vtkOGSTimeStatistics, vtkDataSetAlgorithm);

  void SetStartTime(const char *tstep);
  void SetEndTime(const char *tstep);

  vtkStringArray *GetTimeValues();

protected:
  vtkOGSTimeStatistics();
  ~vtkOGSTimeStatistics();

  int RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int RequestUpdateExtent(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int RequestData(vtkInformation *, vtkInformationVector **,vtkInformationVector *) override;

  void InitializeStatistics(vtkDataSet *input, vtkDataSet *output);
  void AccumulateStatistics(vtkDataSet *input, vtkDataSet *output);
  void FinalizeStatistics(vtkDataSet *input, vtkDataSet *output);

  vtkStringArray *TimeValues;

private:
  vtkOGSTimeStatistics(const vtkOGSTimeStatistics&) = delete;
  void operator=(const vtkOGSTimeStatistics&) = delete;

  int current_time, start_time, end_time;
  int abort;

};

#endif
