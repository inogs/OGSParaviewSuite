/*=========================================================================

  Program:   OGSPlotViews
  Module:    OGSVerticalProfilePlot.h

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef _OGSVerticalProfilePlot_h
#define _OGSVerticalProfilePlot_h

#include "pqSMProxy.h"
#include "pqView.h"
#include "pqPythonView.h"

class OGSVerticalProfilePlot : public pqPythonView
{
  Q_OBJECT
  typedef pqPythonView Superclass;
public:
  static QString OGSVerticalProfilePlotType() { return "OGSVerticalProfilePlot"; }

  /// constructor takes a bunch of init stuff and must have this signature to
  /// satisfy pqView

  OGSVerticalProfilePlot(const QString& group, 
                         const QString& name,
                         vtkSMProxy* renModule, 
                         pqServer* server, 
                         QObject* parent = NULL
                        );
  ~OGSVerticalProfilePlot();

protected:

private:
  OGSVerticalProfilePlot(const OGSVerticalProfilePlot&); // Not implemented.
  void operator=(const OGSVerticalProfilePlot&); // Not implemented.

};

#endif