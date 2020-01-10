/*=========================================================================

  Program:   OGSPlotViews
  Module:    vtkSMOGSMapPlotProxy.h

  Copyright (c) 2019 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef vtkSMOGSMapPlotProxy_h
#define vtkSMOGSMapPlotProxy_h

#include "vtkPVServerManagerRenderingModule.h" //needed for exports

#include "vtkNew.h" // needed for vtkNew.
#include "vtkSMPythonViewProxy.h"

class vtkImageData;
class vtkRenderer;
class vtkSMProxy;
class vtkSMViewProxyInteractorHelper;

class  vtkSMOGSMapPlotProxy : public vtkSMPythonViewProxy
{
public:
  static vtkSMOGSMapPlotProxy* New();
  vtkTypeMacro(vtkSMOGSMapPlotProxy, vtkSMPythonViewProxy);

protected:
  vtkSMOGSMapPlotProxy();
  ~vtkSMOGSMapPlotProxy() override;

  /**
   * Subclasses should override this method to do the actual image capture.
   */
  vtkImageData* CaptureWindowInternal(int magX, int magY) override;

private:
  vtkSMOGSMapPlotProxy(const vtkSMOGSMapPlotProxy&) = delete;
  void operator=(const vtkSMOGSMapPlotProxy&) = delete;
};

#endif
