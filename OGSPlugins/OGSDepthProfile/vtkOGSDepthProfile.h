// -*- c++ -*-
/*=========================================================================

  Program:   OGSDepthProfile
  Module:    vtkOGSDepthProfile.h

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef vtkOGSDepthProfile_h
#define vtkOGSDepthProfile_h

#include "vtkOGSDepthProfile.h"

#include "vtkDataSetAlgorithm.h"
#include "vtkDataSetAttributes.h"

class vtkCharArray;
class vtkPointData;
class vtkAbstractCellLocator;
class vtkStaticCellLocator;

class VTK_EXPORT vtkOGSDepthProfile : public vtkDataSetAlgorithm 
{
public:
  static vtkOGSDepthProfile* New();
  vtkTypeMacro(vtkOGSDepthProfile, vtkDataSetAlgorithm);

  // Description:
  //Specify the data set that will be probed at the input points.
  //The Input gives the geometry (the points and cells) for the output,
  //while the Source is probed (interpolated) to generate the scalars,
  //vectors, etc. for the output points based on the point locations.
  void SetSourceData(vtkDataObject *source);
  vtkDataObject *GetSource();
  void SetSourceConnection(vtkAlgorithmOutput* algOutput);

  // Description:
  //Returns the name of the char array added to the output with values 1 for
  //valid points and 0 for invalid points.
  //Set to "vtkValidPointMask" by default.
  vtkSetStringMacro(ValidPointMaskArrayName)
  vtkGetStringMacro(ValidPointMaskArrayName)

  // Description
  //Set/Get the prototype cell locator to use for probing the source dataset.
  //By default, vtkStaticCellLocator will be used.
   virtual void SetCellLocatorPrototype(vtkAbstractCellLocator*);
   vtkGetObjectMacro(CellLocatorPrototype, vtkAbstractCellLocator);

protected:
  vtkOGSDepthProfile();
  ~vtkOGSDepthProfile() override;

  int RequestInformation(vtkInformation *, vtkInformationVector **, vtkInformationVector *) override;
  int RequestUpdateExtent(vtkInformation *, vtkInformationVector **, vtkInformationVector *) override;
  int RequestData(vtkInformation *, vtkInformationVector **,vtkInformationVector *) override;

  void Interpolate(vtkDataSet *input, vtkDataSet *source, vtkDataSet *output);

  void BuildFieldList(vtkDataSet* source);

  virtual void InitializeForProbing(vtkDataSet *input, vtkDataSet *output);
  virtual void InitializeOutputArrays(vtkPointData *outPD, vtkIdType numPts);

  void DoProbing(vtkDataSet *input, int srcIdx, vtkDataSet *source,vtkDataSet *output);

  vtkCharArray* MaskPoints;

  vtkAbstractCellLocator* CellLocatorPrototype;

  vtkDataSetAttributes::FieldList* CellList;
  vtkDataSetAttributes::FieldList* PointList;

private:
  vtkOGSDepthProfile(const vtkOGSDepthProfile&) = delete;
  void operator=(const vtkOGSDepthProfile&) = delete;

  char *ValidPointMaskArrayName;

  class vtkVectorOfArrays;
  vtkVectorOfArrays* CellArrays;
};

#endif
