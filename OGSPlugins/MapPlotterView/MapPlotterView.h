/*=========================================================================

  Program:   MapPlotterView
  Module:    MapPlotterView.h

  Copyright (c) 2019 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.
=========================================================================*/

#ifndef _MapPlotterView_h
#define _MapPlotterView_h

#include "pqSMProxy.h"
#include "pqView.h"
#include "pqPythonView.h"

class MapPlotterView : public pqPythonView
{
  Q_OBJECT
  typedef pqPythonView Superclass;
public:
  static QString MapPlotterViewType() { return "MapPlotterView"; }

  /// constructor takes a bunch of init stuff and must have this signature to
  /// satisfy pqView

  MapPlotterView(const QString& group, 
                 const QString& name,
                 vtkSMProxy* renModule, 
                 pqServer* server, 
                 QObject* parent = NULL
                );
  ~MapPlotterView();

protected:

private:
  MapPlotterView(const MapPlotterView&); // Not implemented.
  void operator=(const MapPlotterView&); // Not implemented.

};

#endif