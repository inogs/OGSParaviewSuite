/*=========================================================================

  Program:   Projections
  Module:    projections.hpp

  Map projection conversions.

  Copyright (c) 2019 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef PROJECTIONS_H
#define PROJECTIONS_H

namespace PROJ
{

  /*  PROJMERCATOR
    
      Direct Mercator projection
      https://en.wikipedia.org/wiki/Mercator_projection

  */
  void ProjMercator(double lon, double lat, double xy[2]);

  /* PROJINVMERCATOR

    Inverse Mercator projection
    https://en.wikipedia.org/wiki/Mercator_projection

  */
  void ProjInvMercator(double &lon, double &lat, double xy[2]);

  /* PROJCYLINDRICAL

    Cylindrical equirectangular projection
    https://en.wikipedia.org/wiki/Equirectangular_projection

  */
  void ProjCylindrical(double lon, double lat, double xy[2]);

  /* PROJINVCYLINDRICAL

    Cylindrical equirectangular projection
    https://en.wikipedia.org/wiki/Equirectangular_projection

  */
  void ProjInvCylindrical(double lon, double lat, double xy[2]);

}

#endif