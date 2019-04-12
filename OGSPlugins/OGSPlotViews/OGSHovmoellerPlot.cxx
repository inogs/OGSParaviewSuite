/*=========================================================================

  Program:   OGSPlotViews
  Module:    OGSHovmoellerPlot.cxx

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "OGSHovmoellerPlot.h"
#include "vtkSMOGSHovmoellerPlotProxy.h"

//-----------------------------------------------------------------------------
OGSHovmoellerPlot::OGSHovmoellerPlot(const QString& group, 
                                               const QString& name,
                                               vtkSMProxy* renModule, 
                                               pqServer* server, 
                                               QObject* parent)
      : pqPythonView(OGSHovmoellerPlot::OGSHovmoellerPlotType(), group, name, 
         vtkSMOGSHovmoellerPlotProxy::SafeDownCast(renModule), server, parent)

{
}

//-----------------------------------------------------------------------------
OGSHovmoellerPlot::~OGSHovmoellerPlot()
{
}
