// -*- c++ -*-
/*=========================================================================

  Program:   OGSSelectOkuboWeiss
  Module:    vtkOGSSelectOkuboWeiss.h

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef vtkOGSSelectOkuboWeiss_h
#define vtkOGSSelectOkuboWeiss_h

#include "vtkDataSet.h"
#include "vtkThreshold.h"

class vtkDataArraySelection;

class VTK_EXPORT vtkOGSSelectOkuboWeiss : public vtkThreshold {
public:
  static vtkOGSSelectOkuboWeiss* New();
  vtkTypeMacro(vtkOGSSelectOkuboWeiss, vtkThreshold);

  void PrintSelf(ostream& os, vtkIndent indent) VTK_OVERRIDE;

  // Description:
  // Get the name of the mask field to operate
  vtkSetStringMacro(mask_field);

  // Description:
  // The following methods allow selective seleccion of the basins.
  int GetNumberOfOWArrays();
  const char * GetOWArrayName(int index);
  int GetOWArrayIndex(const char* name);
  int GetOWArrayStatus(const char *name);
  void SetOWArrayStatus(const char* name, int status);
  void DisableAllOWArrays();
  void EnableAllOWArrays();

protected:
  vtkOGSSelectOkuboWeiss();
  ~vtkOGSSelectOkuboWeiss() override;

  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *) override;

  vtkDataArraySelection* OWDataArraySelection;

private:
  vtkOGSSelectOkuboWeiss(const vtkOGSSelectOkuboWeiss&) = delete;
  void operator=(const vtkOGSSelectOkuboWeiss&) = delete;

  char *mask_field;
};

#endif
