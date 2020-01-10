/*=========================================================================

  Program:   OGSPlotViews
  Module:    OGSMapPlot.h

  Copyright (c) 2019 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef _OGSMapPlot_h
#define _OGSMapPlot_h

#include "pqSMProxy.h"
#include "pqView.h"
#include "pqPythonView.h"

class OGSMapPlot : public pqPythonView
{
  Q_OBJECT
  typedef pqPythonView Superclass;
public:
  static QString OGSMapPlotType() { return "OGSMapPlot"; }

  /// constructor takes a bunch of init stuff and must have this signature to
  /// satisfy pqView

  OGSMapPlot(const QString& group, 
             const QString& name,
             vtkSMProxy* renModule, 
             pqServer* server, 
             QObject* parent = NULL
            );
  ~OGSMapPlot();

protected:

private:
  OGSMapPlot(const OGSMapPlot&); // Not implemented.
  void operator=(const OGSMapPlot&); // Not implemented.

};

#endif