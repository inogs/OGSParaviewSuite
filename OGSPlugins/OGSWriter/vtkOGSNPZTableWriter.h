// -*- c++ -*-
/*=========================================================================

  Program:   OGSWriter
  Module:    vtkOGSNPZTableWriter.h

  Copyright (c) 2019 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef vtkOGSNPZTableWriter_h
#define vtkOGSNPZTableWriter_h

#include "vtkIOCoreModule.h" // For export macro
#include "vtkWriter.h"

//----------------------------------------------------------------------------

class vtkOGSNPZTableWriter : public vtkWriter {
public:
  static vtkOGSNPZTableWriter *New();
  vtkTypeMacro(vtkOGSNPZTableWriter, vtkWriter);

  //Get / set the filename where data will be stored (when used as a filter).
  vtkSetStringMacro(FileName);
  vtkGetStringMacro(FileName);

  int Write() override; // This is necessary to get Write() wrapped for scripting languages.

protected:
  vtkOGSNPZTableWriter() override;
  ~vtkOGSNPZTableWriter() override;

  int FillInputPortInformation(int, vtkInformation*) override;
  void WriteData() override;

  char *FileName;

private:
  vtkOGSNPZTableWriter(const vtkOGSNPZTableWriter&) = delete;
  void operator=(const vtkOGSNPZTableWriter&) = delete;
};

#endif
