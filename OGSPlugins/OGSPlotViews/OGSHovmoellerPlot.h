/*=========================================================================

  Program:   OGSPlotViews
  Module:    OGSHovmoellerPlot.h

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef _OGSHovmoellerPlot_h
#define _OGSHovmoellerPlot_h

#include "pqSMProxy.h"
#include "pqView.h"
#include "pqPythonView.h"

class OGSHovmoellerPlot : public pqPythonView
{
  Q_OBJECT
  typedef pqPythonView Superclass;
public:
  static QString OGSHovmoellerPlotType() { return "OGSHovmoellerPlot"; }

  /// constructor takes a bunch of init stuff and must have this signature to
  /// satisfy pqView

  OGSHovmoellerPlot(const QString& group, 
                    const QString& name,
                    vtkSMProxy* renModule, 
                    pqServer* server, 
                    QObject* parent = NULL
                   );
  ~OGSHovmoellerPlot();

protected:

private:
  OGSHovmoellerPlot(const OGSHovmoellerPlot&); // Not implemented.
  void operator=(const OGSHovmoellerPlot&); // Not implemented.

};

#endif