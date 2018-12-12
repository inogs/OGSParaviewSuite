/*=========================================================================

  Program:   OGSPlotViews
  Module:    OGSSpaghettiPlot.h

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef _OGSSpaghettiPlot_h
#define _OGSSpaghettiPlot_h

#include "pqSMProxy.h"
#include "pqView.h"
#include "pqPythonView.h"

class OGSSpaghettiPlot : public pqPythonView
{
  Q_OBJECT
  typedef pqPythonView Superclass;
public:
  static QString OGSSpaghettiPlotType() { return "OGSSpaghettiPlot"; }

  /// constructor takes a bunch of init stuff and must have this signature to
  /// satisfy pqView

  OGSSpaghettiPlot(const QString& group, 
                    const QString& name,
                    vtkSMProxy* renModule, 
                    pqServer* server, 
                    QObject* parent = NULL
                   );
  ~OGSSpaghettiPlot();

protected:

private:
  OGSSpaghettiPlot(const OGSSpaghettiPlot&); // Not implemented.
  void operator=(const OGSSpaghettiPlot&); // Not implemented.

};

#endif