/*=========================================================================

  Program:   MapPlotterView
  Module:    MapPlotterView.cxx

  Copyright (c) 2019 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.
=========================================================================*/

#include "MapPlotterView.h"
#include "vtkSMMapPlotterViewProxy.h"

//-----------------------------------------------------------------------------
MapPlotterView::MapPlotterView(const QString& group, 
                               const QString& name,
                               vtkSMProxy* renModule, 
                               pqServer* server, 
                               QObject* parent)
      : pqPythonView(MapPlotterView::MapPlotterViewType(), group, name, 
         vtkSMMapPlotterViewProxy::SafeDownCast(renModule), server, parent)

{
}

//-----------------------------------------------------------------------------
MapPlotterView::~MapPlotterView()
{
}
