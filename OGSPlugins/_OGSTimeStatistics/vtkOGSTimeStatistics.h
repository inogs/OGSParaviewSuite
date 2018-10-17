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

#include "vtkDataArraySelection.h"
#include "vtkDataSetAlgorithm.h"

class VTK_EXPORT vtkOGSTimeStatistics : public vtkDataSetAlgorithm
{
public:
  static vtkOGSTimeStatistics* New();
  vtkTypeMacro(vtkOGSTimeStatistics, vtkDataSetAlgorithm);


protected:
  vtkOGSTimeStatistics();
  ~vtkOGSTimeStatistics();

  int RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int RequestData(vtkInformation *, vtkInformationVector **,vtkInformationVector *) override;

private:
  vtkOGSTimeStatistics(const vtkOGSTimeStatistics&) = delete;
  void operator=(const vtkOGSTimeStatistics&) = delete;

};

#endif
