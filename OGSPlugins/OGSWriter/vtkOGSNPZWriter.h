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

#include <string>

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

  // Set/Get timestep strings
  vtkGetMacro(timeseries, int);
  vtkSetMacro(timeseries, int);
  void SetStartEnd(const int val1, const int val2);

  int Write() override; // This is necessary to get Write() wrapped for scripting languages.
  virtual int ProcessRequest(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;

protected:
  vtkOGSNPZWriter() override;
  ~vtkOGSNPZWriter() override;

  int FillInputPortInformation(int, vtkInformation*) override;
  int RequestData(vtkInformation* , vtkInformationVector** , vtkInformationVector* ) override;
  void WriteData() override;

  char *FileName, *FileName2, *varname;
  double dfact;
  int singlevar, timeseries, append;
  
  std::string path, fnamenoext, ext;

private:
  vtkOGSNPZWriter(const vtkOGSNPZWriter&) = delete;
  void operator=(const vtkOGSNPZWriter&) = delete;

  int ii_start, ii_end, ii_cur;
};

#endif
