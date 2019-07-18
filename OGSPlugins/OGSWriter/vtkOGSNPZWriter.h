// -*- c++ -*-
/*=========================================================================

  Program:   OGSWriter
  Module:    vtkOGSNPZWriter.h

  Copyright (c) 2019 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef vtkOGSNPZWriter_h
#define vtkOGSNPZWriter_h

#include "vtkIOCoreModule.h" // For export macro
#include "vtkWriter.h"

//----------------------------------------------------------------------------

class vtkOGSNPZWriter : public vtkWriter {
public:
  static vtkOGSNPZWriter *New();
  vtkTypeMacro(vtkOGSNPZWriter, vtkWriter);

  //Get / set the filename where data will be stored (when used as a filter).
  vtkSetStringMacro(FileName);
  vtkGetStringMacro(FileName);

  //Get / set the variable name.
  vtkSetStringMacro(varname);
  vtkGetStringMacro(varname);

  // Let the user select a multiplier factor for the depth
  vtkGetMacro(dfact, double);
  vtkSetMacro(dfact, double);

  // Append variables to the file
  vtkGetMacro(append, int);
  vtkSetMacro(append, int);

  // Using just one variable or all
  vtkGetMacro(singlevar, int);
  vtkSetMacro(singlevar, int);

  int Write() override; // This is necessary to get Write() wrapped for scripting languages.

protected:
  vtkOGSNPZWriter();
  ~vtkOGSNPZWriter() override;

  int FillInputPortInformation(int, vtkInformation*) override;
  void WriteData() override;

  char *FileName, *varname;
  double dfact;
  int singlevar, append;

private:
  vtkOGSNPZWriter(const vtkOGSNPZWriter&) = delete;
  void operator=(const vtkOGSNPZWriter&) = delete;
};

#endif
