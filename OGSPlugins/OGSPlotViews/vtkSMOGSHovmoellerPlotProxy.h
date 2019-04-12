/*=========================================================================

  Program:   OGSPlotViews
  Module:    vtkSMOGSHovmoellerPlotProxy.h

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef vtkSMOGSHovmoellerPlotProxy_h
#define vtkSMOGSHovmoellerPlotProxy_h

#include "vtkPVServerManagerRenderingModule.h" //needed for exports

#include "vtkNew.h" // needed for vtkNew.
#include "vtkSMPythonViewProxy.h"

class vtkImageData;
class vtkRenderer;
class vtkSMProxy;
class vtkSMViewProxyInteractorHelper;

class  vtkSMOGSHovmoellerPlotProxy : public vtkSMPythonViewProxy
{
public:
  static vtkSMOGSHovmoellerPlotProxy* New();
  vtkTypeMacro(vtkSMOGSHovmoellerPlotProxy, vtkSMPythonViewProxy);

protected:
  vtkSMOGSHovmoellerPlotProxy();
  ~vtkSMOGSHovmoellerPlotProxy() override;

  /**
   * Subclasses should override this method to do the actual image capture.
   */
  vtkImageData* CaptureWindowInternal(int magX, int magY) override;

private:
  vtkSMOGSHovmoellerPlotProxy(const vtkSMOGSHovmoellerPlotProxy&) = delete;
  void operator=(const vtkSMOGSHovmoellerPlotProxy&) = delete;
};

#endif
