/*=========================================================================

  Program:   OGSPlotViews
  Module:    vtkSMOGSSpaghettiPlotProxy.h

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef vtkSMOGSSpaghettiPlotProxy_h
#define vtkSMOGSSpaghettiPlotProxy_h

#include "vtkPVServerManagerRenderingModule.h" //needed for exports

#include "vtkNew.h" // needed for vtkNew.
#include "vtkSMPythonViewProxy.h"

class vtkImageData;
class vtkRenderer;
class vtkSMProxy;
class vtkSMViewProxyInteractorHelper;

class  vtkSMOGSSpaghettiPlotProxy : public vtkSMPythonViewProxy
{
public:
  static vtkSMOGSSpaghettiPlotProxy* New();
  vtkTypeMacro(vtkSMOGSSpaghettiPlotProxy, vtkSMPythonViewProxy);

protected:
  vtkSMOGSSpaghettiPlotProxy();
  ~vtkSMOGSSpaghettiPlotProxy() override;

  /**
   * Subclasses should override this method to do the actual image capture.
   */
  vtkImageData* CaptureWindowInternal(int magX, int magY) override;

private:
  vtkSMOGSSpaghettiPlotProxy(const vtkSMOGSSpaghettiPlotProxy&) = delete;
  void operator=(const vtkSMOGSSpaghettiPlotProxy&) = delete;
};

#endif
