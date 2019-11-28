/*=========================================================================

  Program:   MapPlotterView
  Module:    vtkMapPlotterView.h

  Copyright (c) 2019 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.
=========================================================================*/

#ifndef vtkMapPlotterView_h
#define vtkMapPlotterView_h

#include "vtkPythonView.h"

#include <map>

class vtkMapPlotterView : public vtkPythonView
{
public:
  static vtkMapPlotterView* New();
  vtkTypeMacro(vtkMapPlotterView, vtkPythonView);

  /*
    Get/Set the Python script.
  */
  vtkSetStringMacro(Script);
  vtkGetStringMacro(Script);

  /*
    Get/Set the variables.
  */
  vtkSetStringMacro(Variables);
  vtkGetStringMacro(Variables);

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
  vtkMapPlotterView();
  ~vtkMapPlotterView() override;

private:
  vtkMapPlotterView(const vtkMapPlotterView&) =delete;
  void operator=(const vtkMapPlotterView&) =delete;

  char *Script, *Variables;
  std::map<std::string, std::string> mapParam;
};

#endif