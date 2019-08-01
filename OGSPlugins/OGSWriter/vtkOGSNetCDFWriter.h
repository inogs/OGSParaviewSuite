// -*- c++ -*-
/*=========================================================================

  Program:   OGSWriter
  Module:    vtkOGSNetCDFWriter.h

  Copyright (c) 2019 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef vtkOGSNetCDFWriter_h
#define vtkOGSNetCDFWriter_h

#include "vtkIOCoreModule.h" // For export macro
#include "vtkWriter.h"

class vtkAbstractArray;
class vtkRectilinearGrid;
class vtkUnstructuredGrid;
class vtkPolyData;

//----------------------------------------------------------------------------

class vtkOGSNetCDFWriter : public vtkWriter {
public:
  static vtkOGSNetCDFWriter *New();
  vtkTypeMacro(vtkOGSNetCDFWriter, vtkWriter);

  //Get / set the filename where data will be stored (when used as a filter).
  vtkSetStringMacro(FileName);
  vtkGetStringMacro(FileName);

  //Get / set the variable name.
  vtkSetStringMacro(varname);
  vtkGetStringMacro(varname);

  //Get / set the OGS file to load.
  vtkSetStringMacro(OGSFile);
  vtkGetStringMacro(OGSFile);

  // Let the user select a multiplier factor for the depth
  vtkGetMacro(dfact, double);
  vtkSetMacro(dfact, double);

  // Using just one variable or all
  vtkGetMacro(singlevar, int);
  vtkSetMacro(singlevar, int);

  // Using just one variable or all
  vtkGetMacro(SaveAll, int);
  vtkSetMacro(SaveAll, int);

  int Write() override; // This is necessary to get Write() wrapped for scripting languages.

protected:
  vtkOGSNetCDFWriter();
  ~vtkOGSNetCDFWriter() override;

  int FillInputPortInformation(int, vtkInformation*) override;
  void WriteData() override;

  char *FileName, *OGSFile, *varname;
  double dfact;
  int singlevar,SaveAll,projId;

private:
  vtkOGSNetCDFWriter(const vtkOGSNetCDFWriter&) = delete;
  void operator=(const vtkOGSNetCDFWriter&) = delete;

  void writeNetCDFRectilinearGrid1Var(vtkRectilinearGrid *, vtkAbstractArray *, bool);
  void writeNetCDFRectilinearGridnVar(vtkRectilinearGrid *, bool);
  void writeNetCDFUnstructuredGrid1Var(vtkUnstructuredGrid *, vtkAbstractArray *, bool);
  void writeNetCDFUnstructuredGridnVar(vtkUnstructuredGrid *, bool);
  void writeNetCDFPolyData1Var(vtkPolyData *, vtkAbstractArray *, bool);
  void writeNetCDFPolyDatanVar(vtkPolyData *, bool);
};

#endif
