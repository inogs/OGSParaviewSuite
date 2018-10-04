// -*- c++ -*-
/*=========================================================================

  Program:   OGSSpatialStats
  Module:    vtkOGSSpatialStatsFromFile.h

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef vtkOGSSpatialStatsFromFile_h
#define vtkOGSSpatialStatsFromFile_h

#include "vtkDataSet.h"
#include "vtkRectilinearGridAlgorithm.h"

class VTK_EXPORT vtkOGSSpatialStatsFromFile : public vtkRectilinearGridAlgorithm
{
public:
  static vtkOGSSpatialStatsFromFile* New();
  vtkTypeMacro(vtkOGSSpatialStatsFromFile, vtkRectilinearGridAlgorithm);

protected:
  vtkOGSSpatialStatsFromFile();
  ~vtkOGSSpatialStatsFromFile() override;

  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *) override;

private:
  vtkOGSSpatialStatsFromFile(const vtkOGSSpatialStatsFromFile&) = delete;
  void operator=(const vtkOGSSpatialStatsFromFile&) = delete;
};

#endif
