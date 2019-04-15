/*=========================================================================

  Program:   OGSPlotViews
  Module:    vtkSMOGSVerticalProfilePlotProxy.h

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef vtkSMOGSVerticalProfilePlotProxy_h
#define vtkSMOGSVerticalProfilePlotProxy_h

#include "vtkPVServerManagerRenderingModule.h" //needed for exports

#include "vtkNew.h" // needed for vtkNew.
#include "vtkSMPythonViewProxy.h"

class vtkImageData;
class vtkRenderer;
class vtkSMProxy;
class vtkSMViewProxyInteractorHelper;

class  vtkSMOGSVerticalProfilePlotProxy : public vtkSMPythonViewProxy
{
public:
  static vtkSMOGSVerticalProfilePlotProxy* New();
  vtkTypeMacro(vtkSMOGSVerticalProfilePlotProxy, vtkSMPythonViewProxy);

protected:
  vtkSMOGSVerticalProfilePlotProxy();
  ~vtkSMOGSVerticalProfilePlotProxy() override;

  /**
   * Subclasses should override this method to do the actual image capture.
   */
  vtkImageData* CaptureWindowInternal(int magX, int magY) override;

private:
  vtkSMOGSVerticalProfilePlotProxy(const vtkSMOGSVerticalProfilePlotProxy&) = delete;
  void operator=(const vtkSMOGSVerticalProfilePlotProxy&) = delete;
};

#endif
