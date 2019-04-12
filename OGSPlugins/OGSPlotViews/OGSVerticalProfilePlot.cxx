/*=========================================================================

  Program:   OGSPlotViews
  Module:    OGSVerticalProfilePlot.cxx

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "OGSVerticalProfilePlot.h"
#include "vtkSMOGSVerticalProfilePlotProxy.h"

//-----------------------------------------------------------------------------
OGSVerticalProfilePlot::OGSVerticalProfilePlot(const QString& group, 
                                               const QString& name,
                                               vtkSMProxy* renModule, 
                                               pqServer* server, 
                                               QObject* parent)
      : pqPythonView(OGSVerticalProfilePlot::OGSVerticalProfilePlotType(), group, name, 
         vtkSMOGSVerticalProfilePlotProxy::SafeDownCast(renModule), server, parent)

{
}

//-----------------------------------------------------------------------------
OGSVerticalProfilePlot::~OGSVerticalProfilePlot()
{
}
