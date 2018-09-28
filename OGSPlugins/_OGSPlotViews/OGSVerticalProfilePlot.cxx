/*=========================================================================

  Program:   Visualization Toolkit
  Module:    OGSVerticalProfilePlot.cxx

  Copyright (c) 
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "OGSVerticalProfilePlot.h"
#include "vtkSMPythonViewProxy.h"

//-----------------------------------------------------------------------------
OGSVerticalProfilePlot::OGSVerticalProfilePlot(const QString& group, 
                                               const QString& name,
                                               vtkSMProxy* renModule, 
                                               pqServer* server, 
                                               QObject* parent)
      : pqPythonView(OGSVerticalProfilePlot::OGSVerticalProfilePlotType(), group, name, 
         vtkSMPythonViewProxy::SafeDownCast(renModule), server, parent)

{
}

//-----------------------------------------------------------------------------
OGSVerticalProfilePlot::~OGSVerticalProfilePlot()
{
}
