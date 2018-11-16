/*=========================================================================

  Program:   OGSPlotViews
  Module:    vtkOGSHovmoellerPlot.h

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef vtkOGSHovmoellerPlot_h
#define vtkOGSHovmoellerPlot_h

#include "vtkDataSet.h"
#include "vtkPythonView.h"

#include <map>

class vtkOGSHovmoellerPlot : public vtkPythonView
{
public:
  static vtkOGSHovmoellerPlot* New();
  vtkTypeMacro(vtkOGSHovmoellerPlot, vtkPythonView);

  /*
    Get/Set the Python script.
  */
  vtkSetStringMacro(Script);
  vtkGetStringMacro(Script);

  /*
  	Set a name-value parameter that will be available to the script
  	when it is run
  */
  void SetParameterInternal(const char* name, const char* value);
  void SetParameter(const char* name, const char* value);
  void SetParameter(const char* name, int value);
  void SetParameter(const char* name, double value);
  void SetParameter(const char* name, double value1, double value2);
  void SetParameter(const char* name, double value1, double value2, double value3);

  /*
    Overrides the base class method to request an addition pass that moves data from the
    server to the client.
  */
  void Update() VTK_OVERRIDE;

protected:
  vtkOGSHovmoellerPlot();
  ~vtkOGSHovmoellerPlot() override;

private:
  vtkOGSHovmoellerPlot(const vtkOGSHovmoellerPlot&) =delete;
  void operator=(const vtkOGSHovmoellerPlot&) =delete;

  char *Script;
  std::map<std::string, std::string> mapParam;
};

#endif