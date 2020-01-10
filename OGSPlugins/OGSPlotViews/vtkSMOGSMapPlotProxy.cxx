/*=========================================================================

  Program:   OGSPlotViews
  Module:    vtkSMOGSMapPlotProxy.cxx

  Copyright (c) 2019 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkSMOGSMapPlotProxy.h"

#include "vtkClientServerStream.h"
#include "vtkDataArray.h"
#include "vtkMultiProcessController.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkProcessModule.h"
#include "vtkOGSMapPlot.h"
#include "vtkRenderWindow.h"
#include "vtkRenderer.h"
#include "vtkSMViewProxyInteractorHelper.h"

vtkStandardNewMacro(vtkSMOGSMapPlotProxy);

//----------------------------------------------------------------------------
vtkSMOGSMapPlotProxy::vtkSMOGSMapPlotProxy() {}

//----------------------------------------------------------------------------
vtkSMOGSMapPlotProxy::~vtkSMOGSMapPlotProxy() {}

//----------------------------------------------------------------------------
vtkImageData* vtkSMOGSMapPlotProxy::CaptureWindowInternal(int magX, int magY) {

  vtkOGSMapPlot* pv = vtkOGSMapPlot::SafeDownCast(
    this->GetClientSideObject());
  if (pv) pv->SetMagnification(magX, magY);

  vtkImageData* image = this->Superclass::CaptureWindowInternal(magX, magY);
  if (pv) pv->SetMagnification(1, 1);

  return image;
}
