/*=========================================================================

  Program:   OGSUtils
  Module:    vtkOGSPointSource.h

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef vtkOGSPointSource_h
#define vtkOGSPointSource_h

#include "vtkFiltersSourcesModule.h" // For export macro
#include "vtkPolyDataAlgorithm.h"

#define VTK_POINT_UNIFORM   1
#define VTK_POINT_SHELL     0

class vtkRandomSequence;

class vtkOGSPointSource : public vtkPolyDataAlgorithm {
public:
  static vtkOGSPointSource *New();
  vtkTypeMacro(vtkOGSPointSource,vtkPolyDataAlgorithm);

  // Description:
  // Set the number of points to generate.
  //vtkSetClampMacro(NumberOfPoints,vtkIdType,1,VTK_ID_MAX);
  //vtkGetMacro(NumberOfPoints,vtkIdType);

  // Description:
  // Set the center of the point cloud.
  vtkSetVector3Macro(Center,double);
  vtkGetVectorMacro(Center,double,3);

  // Description:
  // Set the radius of the point cloud.  If you are
  // generating a Gaussian distribution, then this is
  // the standard deviation for each of x, y, and z.
  vtkSetClampMacro(Radius,double,0.0,VTK_DOUBLE_MAX);
  vtkGetMacro(Radius,double);

  // Description:
  // Lets the user select a multiplier factor for the depth
  vtkGetMacro(DepthScale, double);
  vtkSetMacro(DepthScale, double);

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
  // Lets the user input a depth
  // instead of a point in the mesh
  void SetDepth(double);
  void GetDepth(double&);

  // Description:
  // Specify the distribution to use.  The default is a
  // uniform distribution.  The shell distribution produces
  // random points on the surface of the sphere, none in the interior.
  vtkSetMacro(Distribution,int);
  void SetDistributionToUniform() { this->SetDistribution(VTK_POINT_UNIFORM); };
  void SetDistributionToShell()   { this->SetDistribution(VTK_POINT_SHELL);   };
  vtkGetMacro(Distribution,int);

  // Description:
  // Set/get the desired precision for the output points.
  // vtkAlgorithm::SINGLE_PRECISION - Output single-precision floating point.
  // vtkAlgorithm::DOUBLE_PRECISION - Output double-precision floating point.
  vtkSetMacro(OutputPointsPrecision,int);
  vtkGetMacro(OutputPointsPrecision,int);

  // Description:
  // Set/Get a random sequence generator.
  // By default, the generator in vtkMath is used to maintain backwards
  // compatibility.
  virtual void SetRandomSequence(vtkRandomSequence *randomSequence);
  vtkGetObjectMacro(RandomSequence,vtkRandomSequence);

protected:
  vtkOGSPointSource(vtkIdType numPts=1);
  ~vtkOGSPointSource() override;

  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *) override;

  double Random();

  vtkIdType NumberOfPoints;
  double Center[3], Radius, DepthScale;
  int Distribution, OutputPointsPrecision, Projection;
  vtkRandomSequence* RandomSequence;

private:
  vtkOGSPointSource(const vtkOGSPointSource&) = delete;
  void operator=(const vtkOGSPointSource&) = delete;
};

#endif
