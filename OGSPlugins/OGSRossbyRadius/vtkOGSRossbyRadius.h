/*=========================================================================

  Program:   OGSRossbyRadius
  Module:    vtkOGSRossbyRadius.h

  Copyright (c) 2020 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef vtkOGSRossbyRadius_h
#define vtkOGSRossbyRadius_h

#include "vtkDataSet.h"
#include "vtkRectilinearGridAlgorithm.h"

#include "vtkPVConfig.h" // For PARAVIEW_USE_MPI

#include <vector>

#include "V3.h"
#include "field.h"

#ifdef PARAVIEW_USE_MPI
  class vtkMultiProcessController;
#endif

class VTK_EXPORT vtkOGSRossbyRadius : public vtkRectilinearGridAlgorithm
{
public:
  static vtkOGSRossbyRadius* New();
  vtkTypeMacro(vtkOGSRossbyRadius, vtkRectilinearGridAlgorithm);
  
  // Description:
  // Get the gravity acceleration
  vtkSetMacro(g, double);
  vtkGetMacro(g, double);

  // Description:
  // Get the gravity acceleration
  vtkSetMacro(f_cor_ct, double);
  vtkGetMacro(f_cor_ct, double);
 
  // Description:
  // Decide to use the volume for statistics
  vtkGetMacro(useCtFcor, int);
  vtkSetMacro(useCtFcor, int);
  vtkBooleanMacro(useCtFcor, int);

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
  vtkOGSRossbyRadius();
  ~vtkOGSRossbyRadius() override;

  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *) override;
  int RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;

  #ifdef PARAVIEW_USE_MPI
    vtkMultiProcessController* Controller;
  #endif

private:
  vtkOGSRossbyRadius(const vtkOGSRossbyRadius&) = delete;
  void operator=(const vtkOGSRossbyRadius&) = delete;

  // Depth levels and number of depth levels
  int procId, nProcs;
  double g, f_cor_ct, Omega, epsi;
  std::vector<double> zcoords;
  char *mask_field;
  bool useCtFcor;

  // Auxiliar variables worth conserving
  // they are read once in the RequestInformation and used in the
  // RequestData
  bool iscelld;              // Whether we have cell or point data 
  v3::V3v xyz;               // Stores cell/point coordinates
  field::Field<int> cId2zId; // Cell to depth level connectivity

  bool isReqInfo;            // Set true when request information
};

#endif
