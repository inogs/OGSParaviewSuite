// -*- c++ -*-
/*=========================================================================

  Program:   OGSHovmoeller
  Module:    vtkOGSHovmoeller.h

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef vtkOGSHovmoeller_h
#define vtkOGSHovmoeller_h

#include "vtkDataSetAlgorithm.h"
#include "vtkTableAlgorithm.h"
#include "vtkDataSetAttributes.h"

class vtkTable;
class vtkStringArray;
class vtkAbstractCellLocator;

class VTK_EXPORT vtkOGSHovmoeller : public vtkDataSetAlgorithm 
{
public:
  static vtkOGSHovmoeller* New();
  vtkTypeMacro(vtkOGSHovmoeller, vtkDataSetAlgorithm);

  // Description:
  //Specify the data set that will be probed at the input points.
  //The Input gives the geometry (the points and cells) for the output,
  //while the Source is probed (interpolated) to generate the scalars,
  //vectors, etc. for the output points based on the point locations.
  void SetSourceData(vtkDataObject *source);
  vtkDataObject *GetSource();
  void SetSourceConnection(vtkAlgorithmOutput* algOutput);

  // Description
  //Set/Get the prototype cell locator to use for probing the source dataset.
  //By default, vtkStaticCellLocator will be used.
   virtual void SetCellLocatorPrototype(vtkAbstractCellLocator*);
   vtkGetObjectMacro(CellLocatorPrototype, vtkAbstractCellLocator);

  void SetStartTime(const char *tstep);
  void SetEndTime(const char *tstep);

  vtkStringArray *GetTimeValues();

protected:
  vtkOGSHovmoeller();
  ~vtkOGSHovmoeller() override;

  int RequestInformation(vtkInformation *, vtkInformationVector **, vtkInformationVector *) override;
  int RequestUpdateExtent(vtkInformation *, vtkInformationVector **, vtkInformationVector *) override;
  int RequestData(vtkInformation *, vtkInformationVector **,vtkInformationVector *) override;

//  void Initialize(vtkDataSet* input,vtkDataSet* source, vtkTable* output);
  void Interpolate(vtkDataSet *input, vtkDataSet *source, vtkTable *output);

  vtkStringArray *TimeValues;

private:
  vtkOGSHovmoeller(const vtkOGSHovmoeller&) = delete;
  void operator=(const vtkOGSHovmoeller&) = delete;

  vtkAbstractCellLocator* CellLocatorPrototype;
  
  vtkDataSetAttributes::FieldList* CellList;
  vtkDataSetAttributes::FieldList* PointList;

  int current_time, start_time, end_time;
  int abort;
};

#endif
