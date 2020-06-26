/*=========================================================================

  Program:   OGSMixingLayerDepth
  Module:    vtkOGSMixingLayerDepth.h

  Copyright (c) 2020 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef vtkOGSMixingLayerDepth_h
#define vtkOGSMixingLayerDepth_h

#include "vtkDataSet.h"
#include "vtkRectilinearGridAlgorithm.h"

#include "vtkPVConfig.h" // For PARAVIEW_USE_MPI

#include <vector>

#include "V3.h"
#include "field.h"

#ifdef PARAVIEW_USE_MPI
  class vtkMultiProcessController;
#endif

class VTK_EXPORT vtkOGSMixingLayerDepth : public vtkRectilinearGridAlgorithm
{
public:
  static vtkOGSMixingLayerDepth* New();
  vtkTypeMacro(vtkOGSMixingLayerDepth, vtkRectilinearGridAlgorithm);
  
  // Description:
  // Get the temperature difference
  vtkSetMacro(dT, double);
  vtkGetMacro(dT, double);
 
  // Description:
  // Decide to use the volume for statistics
  vtkGetMacro(useDensity, int);
  vtkSetMacro(useDensity, int);
  vtkBooleanMacro(useDensity, int);

  // Description:
  // Get the density difference
  vtkSetMacro(drho, double);
  vtkGetMacro(drho, double);
  
  // Description:
  // Get the reference depth
  vtkSetMacro(zref, double);
  vtkGetMacro(zref, double);

  // Description:
  // Get the epsilon
  vtkSetMacro(epsi, double);
  vtkGetMacro(epsi, double);

  // Description:
  // Get the name of the mask field to operate
  vtkSetStringMacro(mask_field);

  #ifdef PARAVIEW_USE_MPI
    // Description:
    // Set the controller use in compositing (set to
    // the global controller by default)
    // If not using the default, this must be called before any
    // other methods.
    virtual void SetController(vtkMultiProcessController* controller);
  #endif

protected:
  vtkOGSMixingLayerDepth();
  ~vtkOGSMixingLayerDepth() override;

  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *) override;
  int RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;

  #ifdef PARAVIEW_USE_MPI
    vtkMultiProcessController* Controller;
  #endif

private:
  vtkOGSMixingLayerDepth(const vtkOGSMixingLayerDepth&) = delete;
  void operator=(const vtkOGSMixingLayerDepth&) = delete;

  // Depth levels and number of depth levels
  int method, procId, nProcs;
  double dT, drho, zref, epsi;
  std::vector<double> zcoords;
  bool useDensity;
  char *mask_field;

  // Auxiliar variables worth conserving
  // they are read once in the RequestInformation and used in the
  // RequestData
  bool iscelld;              // Whether we have cell or point data 
  v3::V3v xyz;               // Stores cell/point coordinates
  field::Field<int> cId2zId; // Cell to depth level connectivity

  bool isReqInfo;            // Set true when request information
};

#endif
