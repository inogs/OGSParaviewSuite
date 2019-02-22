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

#include <vector>

#include "../_utils/V3.h"
#include "../_utils/field.h"

class VTK_EXPORT vtkOGSSpatialStats : public vtkDataSetAlgorithm
{
public:
  static vtkOGSSpatialStats* New();
  vtkTypeMacro(vtkOGSSpatialStats, vtkDataSetAlgorithm);

  // Description:
  // Get the epsilon
  vtkSetMacro(epsi, double);
  vtkGetMacro(epsi, double);

  // Description:
  // Number of user inputed depths
  void SetNumberOfDepthLevels(int n);
  void SetDepthLevels(int i, double value);

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
  int RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;

  void CellStats(vtkDataSet *, vtkDataSet *, double );
  void PointStats(vtkDataSet *, vtkDataSet *, double );

  vtkDataArraySelection* StatDataArraySelection;

private:
  vtkOGSSpatialStats(const vtkOGSSpatialStats&) = delete;
  void operator=(const vtkOGSSpatialStats&) = delete;

  // Tolerance for finding the depth levels
  double epsi;

  // Depth levels and number of depth levels
  int ndepths;
  std::vector<double> zcoords;

  // Auxiliar variables worth conserving
  // they are read once in the RequestInformation and used in the
  // RequestData
  bool iscelld;              // Whether we have cell or point data 
  v3::V3v xyz;               // Stores cell/point coordinates
  field::Field<int> cId2zId; // Cell to depth level connectivity

  bool isReqInfo;            // Set true when request information
};

#endif
