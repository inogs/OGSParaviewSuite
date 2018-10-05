// -*- c++ -*-
/*=========================================================================

  Program:   OGSSpatialStats
  Module:    vtkOGSSpatialStats.h

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef vtkOGSSpatialStats_h
#define vtkOGSSpatialStats_h

#include "vtkDataSet.h"

#include "vtkDataArraySelection.h"
#include "vtkDataSetAlgorithm.h"

class VTK_EXPORT vtkOGSSpatialStats : public vtkDataSetAlgorithm
{
public:
  static vtkOGSSpatialStats* New();
  vtkTypeMacro(vtkOGSSpatialStats, vtkDataSetAlgorithm);

  // Description:
  // The following methods allow selective seleccion of aggregate variables.
  int GetNumberOfStatArrays();
  const char * GetStatArrayName(int index);
  int GetStatArrayIndex(const char* name);
  int GetStatArrayStatus(const char *name);
  void SetStatArrayStatus(const char* name, int status);
  void DisableAllStatArrays();
  void EnableAllStatArrays();

protected:
  vtkOGSSpatialStats();
  ~vtkOGSSpatialStats() override;

  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *) override;

  vtkDataArraySelection* StatDataArraySelection;

private:
  vtkOGSSpatialStats(const vtkOGSSpatialStats&) = delete;
  void operator=(const vtkOGSSpatialStats&) = delete;

};

#endif
