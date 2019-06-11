/*=========================================================================

  Program:   Projections
  Module:    projections.cpp

  Map projection conversions.

  Copyright (c) 2019 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include<cmath>
#include "Projections.hpp"

// Usefil definitions
#define PI             4*std::atan(1.)
#define DEG2RAD(theta) PI*(theta)/180.
#define RAD2DEG(theta) 180.*(theta)/PI

// Some user defined constants to work with the Python modules
// of this program
const double R_earth = 6371.e3; // Earth mean radius (m)
const double lon0    = -5.3;
const double lat0    = 28.;

namespace PROJ
{
	/*	PROJMERCATOR
		
		Direct Mercator projection
		https://en.wikipedia.org/wiki/Mercator_projection

	*/
	void ProjMercator(double lon, double lat, double xy[2]) {

		xy[0] = R_earth*( DEG2RAD(lon-lon0) ); // x coordinate

		double mercN  = std::log( std::tan(DEG2RAD(45. +  lat/2.)) );
		double mercN0 = std::log( std::tan(DEG2RAD(45. + lat0/2.)) );
		xy[1] = R_earth*(mercN - mercN0);      // y coordinate
	}

	/* PROJINVMERCATOR

		Inverse Mercator projection
		https://en.wikipedia.org/wiki/Mercator_projection

	*/
	void ProjInvMercator(double &lon, double &lat, double xy[2]) {

		lon = lon0 + RAD2DEG(xy[0]/R_earth);

		double mercN0 = std::log( std::tan(DEG2RAD(45. + lat0/2.)) );
		lat = RAD2DEG( 2.*std::atan(std::exp(xy[1]/R_earth + mercN0)) - PI/2. );
	}

	/* PROJCYLINDRICAL

		Cylindrical equirectangular projection
		https://en.wikipedia.org/wiki/Equirectangular_projection

	*/
	void ProjCylindrical(double lon, double lat, double xy[2]) {

		xy[0] = R_earth*DEG2RAD(lon);
		xy[1] = R_earth*DEG2RAD(lat);
	}

	/* PROJINVCYLINDRICAL

		Cylindrical equirectangular projection
		https://en.wikipedia.org/wiki/Equirectangular_projection

	*/
	void ProjInvCylindrical(double lon, double lat, double xy[2]) {

		lon = RAD2DEG(xy[0]/R_earth);
		lat = RAD2DEG(xy[1]/R_earth);
	}
}