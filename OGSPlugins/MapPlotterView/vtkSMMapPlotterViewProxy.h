/*=========================================================================

  Program:   MapPlotterView
  Module:    vtkSMMapPlotterViewProxy.h

  Copyright (c) 2019 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.
=========================================================================*/

#ifndef vtkSMMapPlotterViewProxy_h
#define vtkSMMapPlotterViewProxy_h

#include "vtkPVServerManagerRenderingModule.h" //needed for exports

#include "vtkNew.h" // needed for vtkNew.
#include "vtkSMPythonViewProxy.h"

class vtkImageData;
class vtkRenderer;
class vtkSMProxy;
class vtkSMViewProxyInteractorHelper;

class VTKPVSERVERMANAGERRENDERING_EXPORT vtkSMMapPlotterViewProxy : public vtkSMPythonViewProxy
{
public:
  static vtkSMMapPlotterViewProxy* New();
  vtkTypeMacro(vtkSMMapPlotterViewProxy, vtkSMPythonViewProxy);

protected:
  vtkSMMapPlotterViewProxy();
  ~vtkSMMapPlotterViewProxy() override;

  /**
   * Subclasses should override this method to do the actual image capture.
   */
  vtkImageData* CaptureWindowInternal(int magX, int magY) override;

private:
  vtkSMMapPlotterViewProxy(const vtkSMMapPlotterViewProxy&) = delete;
  void operator=(const vtkSMMapPlotterViewProxy&) = delete;
};

#endif
