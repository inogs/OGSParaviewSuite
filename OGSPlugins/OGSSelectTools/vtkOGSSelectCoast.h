// -*- c++ -*-
/*=========================================================================

  Program:   OGSSelectCoast
  Module:    vtkOGSSelectCoast.h

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef vtkOGSSelectCoast_h
#define vtkOGSSelectCoast_h

#include "vtkDataSet.h"
#include "vtkThreshold.h"

class vtkDataArraySelection;

class VTK_EXPORT vtkOGSSelectCoast : public vtkThreshold {
public:
  static vtkOGSSelectCoast* New();
  vtkTypeMacro(vtkOGSSelectCoast, vtkThreshold);

  void PrintSelf(ostream& os, vtkIndent indent) VTK_OVERRIDE;

  // Description:
  // Get the name of the mask field to operate
  vtkSetStringMacro(mask_field);

  // Description:
  // The following methods allow selective seleccion of the basins.
  int GetNumberOfCoastsArrays();
  const char * GetCoastsArrayName(int index);
  int GetCoastsArrayIndex(const char* name);
  int GetCoastsArrayStatus(const char *name);
  void SetCoastsArrayStatus(const char* name, int status);
  void DisableAllCoastsArrays();
  void EnableAllCoastsArrays();

protected:
  vtkOGSSelectCoast();
  ~vtkOGSSelectCoast();

  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *) override;

  vtkDataArraySelection* CoastsDataArraySelection;

private:
  vtkOGSSelectCoast(const vtkOGSSelectCoast&) = delete;
  void operator=(const vtkOGSSelectCoast&) = delete;

  char *mask_field;
};

#endif
