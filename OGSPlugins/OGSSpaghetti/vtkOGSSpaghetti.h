// -*- c++ -*-
/*=========================================================================

  Program:   OGSSpaghetti
  Module:    vtkOGSSpaghetti.h

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef vtkOGSSpaghetti_h
#define vtkOGSSpaghetti_h

#include "vtkDataSetAlgorithm.h"
#include "vtkTableAlgorithm.h"
#include "vtkDataSetAttributes.h"

class vtkTable;
class vtkStringArray;
class vtkAbstractCellLocator;

class VTK_EXPORT vtkOGSSpaghetti : public vtkDataSetAlgorithm 
{
public:
  static vtkOGSSpaghetti* New();
  vtkTypeMacro(vtkOGSSpaghetti, vtkDataSetAlgorithm);

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

  // Description:
  // Get the name of the variable field
  vtkSetStringMacro(field);

  // Description:
  // Lets the user select a multiplier factor for the depth
  vtkGetMacro(DepthScale, double);
  vtkSetMacro(DepthScale, double);

protected:
  vtkOGSSpaghetti();
  ~vtkOGSSpaghetti() override;

  int FillOutputPortInformation(int , vtkInformation*);
  int RequestDataObject(vtkInformation *, vtkInformationVector **, vtkInformationVector *) override;
  int RequestInformation(vtkInformation *, vtkInformationVector **, vtkInformationVector *) override;
  int RequestUpdateExtent(vtkInformation *, vtkInformationVector **, vtkInformationVector *) override;
  int RequestData(vtkInformation *, vtkInformationVector **,vtkInformationVector *) override;

//  void Initialize(vtkDataSet* input,vtkDataSet* source, vtkTable* output);
  void Interpolate(vtkDataSet *input, vtkDataSet *source, vtkTable *output);

  vtkStringArray *TimeValues;

private:
  vtkOGSSpaghetti(const vtkOGSSpaghetti&) = delete;
  void operator=(const vtkOGSSpaghetti&) = delete;

  vtkAbstractCellLocator* CellLocatorPrototype;
  
  vtkDataSetAttributes::FieldList* CellList;
  vtkDataSetAttributes::FieldList* PointList;

  int current_time, start_time, end_time;
  int abort;

  char *field;
  double DepthScale;
};

#endif
