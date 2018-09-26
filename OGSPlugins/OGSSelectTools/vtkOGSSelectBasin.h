// -*- c++ -*-
/*=========================================================================

  Program:   OGSSelectBasin
  Module:    vtkOGSSelectBasin.h

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef vtkOGSSelectBasin_h
#define vtkOGSSelectBasin_h

#include "vtkDataSet.h"
#include "vtkThreshold.h"

class vtkDataArraySelection;

class VTK_EXPORT vtkOGSSelectBasin : public vtkThreshold {
public:
  static vtkOGSSelectBasin* New();
  vtkTypeMacro(vtkOGSSelectBasin, vtkThreshold);

  void PrintSelf(ostream& os, vtkIndent indent) VTK_OVERRIDE;

  // Description:
  // Get the name of the mask field to operate
  vtkSetStringMacro(mask_field);

  // Description:
  // The following methods allow selective seleccion of the basins.
  int GetNumberOfBasinsArrays();
  const char * GetBasinsArrayName(int index);
  int GetBasinsArrayIndex(const char* name);
  int GetBasinsArrayStatus(const char *name);
  void SetBasinsArrayStatus(const char* name, int status);
  void DisableAllBasinsArrays();
  void EnableAllBasinsArrays();

protected:
  vtkOGSSelectBasin();
  ~vtkOGSSelectBasin() override;

  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *) override;

  vtkDataArraySelection* BasinsDataArraySelection;

private:
  vtkOGSSelectBasin(const vtkOGSSelectBasin&) = delete;
  void operator=(const vtkOGSSelectBasin&) = delete;

  char *mask_field;
};

#endif
