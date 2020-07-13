// -*- c++ -*-
/*=========================================================================

  Program:   OGSWriter
  Module:    vtkOGSFieldWriter.h

  Copyright (c) 2019 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef vtkOGSFieldWriter_h
#define vtkOGSFieldWriter_h

#include "vtkIOCoreModule.h" // For export macro
#include "vtkWriter.h"

#include <string>

//----------------------------------------------------------------------------

class vtkOGSFieldWriter : public vtkWriter {
public:
  static vtkOGSFieldWriter *New();
  vtkTypeMacro(vtkOGSFieldWriter, vtkWriter);

  //Get / set the filename where data will be stored (when used as a filter).
  vtkSetStringMacro(FileName);
  vtkGetStringMacro(FileName);

  //Get / set the variable name.
  vtkSetStringMacro(varname);
  vtkGetStringMacro(varname);

  // Let the user select a multiplier factor for the depth
  vtkGetMacro(dfact, double);
  vtkSetMacro(dfact, double);

  // Use Fortran or C/C++ stride
  vtkGetMacro(fstride, int);
  vtkSetMacro(fstride, int);

  // Set/Get timestep strings
  vtkGetMacro(timeseries, int);
  vtkSetMacro(timeseries, int);
  void SetStartEnd(const int val1, const int val2);

  int Write() override; // This is necessary to get Write() wrapped for scripting languages.
  virtual int ProcessRequest(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;

protected:
  vtkOGSFieldWriter() override;
  ~vtkOGSFieldWriter() override;

  int FillInputPortInformation(int, vtkInformation*) override;
  int RequestData(vtkInformation* , vtkInformationVector** , vtkInformationVector* ) override;
  void WriteData() override;

  char *FileName, *varname;
  double dfact;
  int fstride, timeseries;

  std::string path, fnamenoext, ext;

private:
  vtkOGSFieldWriter(const vtkOGSFieldWriter&) = delete;
  void operator=(const vtkOGSFieldWriter&) = delete;

  int ii_start, ii_end, ii_cur;
};

#endif
