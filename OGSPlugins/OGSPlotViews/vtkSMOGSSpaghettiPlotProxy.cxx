/*=========================================================================

  Program:   OGSPlotViews
  Module:    vtkSMOGSSpaghettiPlotProxy.cxx

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkSMOGSSpaghettiPlotProxy.h"

#include "vtkClientServerStream.h"
#include "vtkDataArray.h"
#include "vtkMultiProcessController.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkProcessModule.h"
#include "vtkOGSSpaghettiPlot.h"
#include "vtkRenderWindow.h"
#include "vtkRenderer.h"
#include "vtkSMViewProxyInteractorHelper.h"

vtkStandardNewMacro(vtkSMOGSSpaghettiPlotProxy);

//----------------------------------------------------------------------------
vtkSMOGSSpaghettiPlotProxy::vtkSMOGSSpaghettiPlotProxy() {}

//----------------------------------------------------------------------------
vtkSMOGSSpaghettiPlotProxy::~vtkSMOGSSpaghettiPlotProxy() {}

//----------------------------------------------------------------------------
vtkImageData* vtkSMOGSSpaghettiPlotProxy::CaptureWindowInternal(int magX, int magY) {

  vtkOGSSpaghettiPlot* pv = vtkOGSSpaghettiPlot::SafeDownCast(
    this->GetClientSideObject());
  if (pv) pv->SetMagnification(magX, magY);

  vtkImageData* image = this->Superclass::CaptureWindowInternal(magX, magY);
  if (pv) pv->SetMagnification(1, 1);

  return image;
}
