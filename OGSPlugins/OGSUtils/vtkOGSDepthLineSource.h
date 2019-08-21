/*=========================================================================

  Program:   OGSUtils
  Module:    vtkOGSDepthLineSource.h

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef vtkOGSDepthLineSource_h
#define vtkOGSDepthLineSource_h

#include "vtkFiltersSourcesModule.h" // For export macro
#include "vtkPolyDataAlgorithm.h"

class vtkPoints;

class vtkOGSDepthLineSource : public vtkPolyDataAlgorithm
{
public:
  static vtkOGSDepthLineSource *New();
  vtkTypeMacro(vtkOGSDepthLineSource,vtkPolyDataAlgorithm);

  // Description:
  // Set/Get position of first end point.
  void SetPoint1(double,double,double);
  void SetPoint1(double[3]);
  void SetPoint1(float[3]);
  void GetPoint1(double&,double&,double&);
  void GetPoint1(double[3]);

  // Description:
  // Set/Get position of other end point.
  void SetPoint2(double,double,double);
  void SetPoint2(double[3]);
  void SetPoint2(float[3]);
  void GetPoint2(double&,double&,double&);
  void GetPoint2(double[3]);

  // Description:
  // Lets the user select a multiplier factor for the depth
  vtkGetMacro(DepthScale, double);
  void SetDepthScale(double);
//  vtkSetMacro(DepthScale, double);

  // Description:
  // Lets the user select a projection type
  vtkGetMacro(Projection, int);
  vtkSetMacro(Projection, int);

  // Description:
  // Lets the user input a Longitude and Latitude pair 
  // instead of a point in the mesh
  void SetLonLat(double,double);
  void GetLonLat(double&,double&);

  // Description:
  // Lets the user input a depth range
  void SetDepthRange(double,double);
  void GetDepthRange(double&,double&);

  // Description:
  // Set/Get the list of points defining a broken line
  virtual void SetPoints(vtkPoints*);
  vtkGetObjectMacro(Points,vtkPoints);

  // Description:
  // Divide line into Resolution number of pieces.
  vtkSetClampMacro(Resolution,int,1,VTK_INT_MAX);
  vtkGetMacro(Resolution,int);

  // Description:
  // Set/get the desired precision for the output points.
  // vtkAlgorithm::SINGLE_PRECISION - Output single-precision floating point.
  // vtkAlgorithm::DOUBLE_PRECISION - Output double-precision floating point.
  vtkSetMacro(OutputPointsPrecision,int);
  vtkGetMacro(OutputPointsPrecision,int);

protected:
  vtkOGSDepthLineSource(int res=1);
  ~vtkOGSDepthLineSource() override;

  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *) override;
  int RequestInformation(vtkInformation *, vtkInformationVector **, vtkInformationVector *) override;

  double Point1[3], Point2[3], DepthScale;
  int Resolution, OutputPointsPrecision, Projection;

  
  // The list of points defining a broken line
  // NB: The Point1/Point2 definition of a single line segment is used by default
  vtkPoints* Points;

private:
  vtkOGSDepthLineSource(const vtkOGSDepthLineSource&) = delete;
  void operator=(const vtkOGSDepthLineSource&) = delete;
};

#endif
