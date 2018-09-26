/*=========================================================================

  Program:   ParaView
  Module:    vtkOGSVerticalProfilePlot.h

=========================================================================*/

#ifndef vtkOGSVerticalProfilePlot_h
#define vtkOGSVerticalProfilePlot_h

#include "vtkDataSet.h"
#include "vtkPythonView.h"

class vtkOGSVerticalProfilePlot : public vtkPythonView
{
public:
  static vtkOGSVerticalProfilePlot* New();
  vtkTypeMacro(vtkOGSVerticalProfilePlot, vtkPythonView);

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
  vtkOGSVerticalProfilePlot();
  ~vtkOGSVerticalProfilePlot() override;

private:
  vtkOGSVerticalProfilePlot(const vtkOGSVerticalProfilePlot&) =delete;
  void operator=(const vtkOGSVerticalProfilePlot&) =delete;

  char *Script;
  char Params[512];
};

#endif