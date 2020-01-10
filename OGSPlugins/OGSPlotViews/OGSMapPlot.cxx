/*=========================================================================

  Program:   OGSPlotViews
  Module:    OGSMapPlot.cxx

  Copyright (c) 2019 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "OGSMapPlot.h"
#include "vtkSMOGSMapPlotProxy.h"

//-----------------------------------------------------------------------------
OGSMapPlot::OGSMapPlot(const QString& group, 
                       const QString& name,
                       vtkSMProxy* renModule, 
                       pqServer* server, 
                       QObject* parent)
      : pqPythonView(OGSMapPlot::OGSMapPlotType(), group, name, 
         vtkSMOGSMapPlotProxy::SafeDownCast(renModule), server, parent)

{
}

//-----------------------------------------------------------------------------
OGSMapPlot::~OGSMapPlot()
{
}
