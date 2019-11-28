/*=========================================================================

  Program:   MapPlotterView
  Module:    vtkSMMapPlotterViewProxy.cxx

  Copyright (c) 2019 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.
=========================================================================*/

#include "vtkSMMapPlotterViewProxy.h"

#include "vtkClientServerStream.h"
#include "vtkDataArray.h"
#include "vtkMultiProcessController.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkProcessModule.h"
#include "vtkMapPlotterView.h"
#include "vtkRenderWindow.h"
#include "vtkRenderer.h"
#include "vtkSMViewProxyInteractorHelper.h"

vtkStandardNewMacro(vtkSMMapPlotterViewProxy);

//----------------------------------------------------------------------------
vtkSMMapPlotterViewProxy::vtkSMMapPlotterViewProxy()
{
}

//----------------------------------------------------------------------------
vtkSMMapPlotterViewProxy::~vtkSMMapPlotterViewProxy()
{
}

//----------------------------------------------------------------------------
vtkImageData* vtkSMMapPlotterViewProxy::CaptureWindowInternal(int magX, int magY)
{
  vtkMapPlotterView* pv = vtkMapPlotterView::SafeDownCast(this->GetClientSideObject());
  if (pv)
  {
    pv->SetMagnification(magX, magY);
  }
  vtkImageData* image = this->Superclass::CaptureWindowInternal(magX, magY);
  if (pv)
  {
    pv->SetMagnification(1, 1);
  }
  return image;
}
