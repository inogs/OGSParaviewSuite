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

#include "vtkTableAlgorithm.h"
#include "vtkDataSetAttributes.h"
#include "vtkPVConfig.h" // For PARAVIEW_USE_MPI

#include <vector>
#include <string>
#include "TimeInterval.h"
#include "TimeList.h"

#include "V3.h"
#include "field.h"

class vtkStringArray;
class vtkAbstractCellLocator;

#ifdef PARAVIEW_USE_MPI
class vtkMultiProcessController;
#endif

class VTK_EXPORT vtkOGSHovmoeller : public vtkTableAlgorithm 
{
public:
  static vtkOGSHovmoeller* New();
  vtkTypeMacro(vtkOGSHovmoeller, vtkTableAlgorithm);

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

  void SetStartTI(const char *tstep);
  void SetEndTI(const char *tstep);

  // Description:
  // Get the name of the variable field
  vtkSetStringMacro(field);
  vtkGetStringMacro(field);
  
  // Description:
  // Averaging algorithm
  vtkGetMacro(use_average, int);
  vtkSetMacro(use_average, int);
  vtkBooleanMacro(use_average, int);

  vtkGetMacro(use_files, int);
  vtkSetMacro(use_files, int);
  vtkBooleanMacro(use_files, int);

  // Description:
  // Folder to STATE_PROFILES
  vtkSetStringMacro(FolderName);
  vtkGetStringMacro(FolderName);

  // Description:
  // Statistic to plot
  vtkGetMacro(sId, int);
  vtkSetMacro(sId, int);

  // Description:
  // Statistics per coast
  vtkGetMacro(per_coast, int);
  vtkSetMacro(per_coast, int);
  vtkBooleanMacro(per_coast, int);

  // Description:
  // Get the epsilon
  vtkSetMacro(epsi, double);
  vtkGetMacro(epsi, double);

  // Description:
  // Get the name of the mask fields to operate
  vtkSetStringMacro(bmask_field);
  vtkSetStringMacro(cmask_field);

  vtkStringArray *GetTimeValues();

  #ifdef PARAVIEW_USE_MPI
  vtkGetObjectMacro(Controller, vtkMultiProcessController);
  virtual void SetController(vtkMultiProcessController*);
  #endif

protected:
  vtkOGSHovmoeller();
  ~vtkOGSHovmoeller() override;

  int FillInputPortInformation(int , vtkInformation *) override;
  int RequestInformation(vtkInformation *, vtkInformationVector **, vtkInformationVector *) override;
  int RequestUpdateExtent(vtkInformation *, vtkInformationVector **, vtkInformationVector *) override;
  int RequestData(vtkInformation *, vtkInformationVector **,vtkInformationVector *) override;

  #ifdef PARAVIEW_USE_MPI
  vtkMultiProcessController* Controller;
  #endif

private:
  vtkOGSHovmoeller(const vtkOGSHovmoeller&) = delete;
  void operator=(const vtkOGSHovmoeller&) = delete;

  int PipelineIterationAlgorithm(vtkInformation *, vtkDataSet *, vtkDataSet *, vtkTable *);
  int FileIterationAlgorithm(vtkInformation *, vtkDataSet *, vtkDataSet *, vtkTable *);
  int AveragesIterationAlgorithm(int, vtkInformation *, vtkDataSet *, vtkDataSet *, vtkTable *);

  vtkAbstractCellLocator* CellLocatorPrototype;
  
  Time::TimeInterval TI;       // TimeInterval for the generic requestor
  Time::TimeList TL;           // TimeList containing all the instants

  std::vector<int> instants;   // Instant ID to loop
  std::vector<double> weights; // Weights for the instants

  int procId, nProcs;
  int CurrentTimeIndex, sId, per_coast;

  char *field, *FolderName, *bmask_field, *cmask_field;
  bool TL_computed, isReqInfo, use_average, use_files;

  std::string tstep_st, tstep_ed;

  double epsi, dfact;

  std::vector<double> zcoords;
  v3::V3v xyz;               // Stores cell/point coordinates
  field::Field<int> cId2zId; // Cell to depth level connectivity
};

#endif
