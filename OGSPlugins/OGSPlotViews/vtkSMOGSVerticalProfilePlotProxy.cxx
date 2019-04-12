/*=========================================================================

  Program:   OGSPlotViews
  Module:    vtkSMOGSVerticalProfilePlotProxy.cxx

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkSMOGSVerticalProfilePlotProxy.h"

#include "vtkClientServerStream.h"
#include "vtkDataArray.h"
#include "vtkMultiProcessController.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkProcessModule.h"
#include "vtkOGSVerticalProfilePlot.h"
#include "vtkRenderWindow.h"
#include "vtkRenderer.h"
#include "vtkSMViewProxyInteractorHelper.h"

vtkStandardNewMacro(vtkSMOGSVerticalProfilePlotProxy);

//----------------------------------------------------------------------------
vtkSMOGSVerticalProfilePlotProxy::vtkSMOGSVerticalProfilePlotProxy() {}

//----------------------------------------------------------------------------
vtkSMOGSVerticalProfilePlotProxy::~vtkSMOGSVerticalProfilePlotProxy() {}

//----------------------------------------------------------------------------
vtkImageData* vtkSMOGSVerticalProfilePlotProxy::CaptureWindowInternal(int magX, int magY) {

  vtkOGSVerticalProfilePlot* pv = vtkOGSVerticalProfilePlot::SafeDownCast(
    this->GetClientSideObject());
  if (pv) pv->SetMagnification(magX, magY);

  vtkImageData* image = this->Superclass::CaptureWindowInternal(magX, magY);
  if (pv) pv->SetMagnification(1, 1);

  return image;
}
