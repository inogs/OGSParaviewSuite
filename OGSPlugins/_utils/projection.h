/*=========================================================================

	Module:    Projection

	Map projection conversions using the proj API.

	Copyright (c) 2019 Arnau Miro, OGS
	All rights reserved.

		 This software is distributed WITHOUT ANY WARRANTY; without even
		 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
		 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef PROJECTIONS_H
#define PROJECTIONS_H

#include <string>
#include <map>

#include "proj_api.h"

#define PROJ_ERR 0
#define PROJ_OK  1

namespace PROJ
{
	inline double rad2deg(const double a) { return (180.0/3.14159265359)*a; }
	inline double deg2rad(const double a) { return (3.14159265359/180.0)*a; }

	class Projection {
		/* PROJECTION

			Class interface with the proj library
		*/
		public:
			inline Projection()  {};
			inline ~Projection() { free(); };

			inline void add(const char* key, const char* projstr) { projs[key] = projstr; }

			inline int transform_point(const char *from, const char *to, double xy[], double out[]) { 
				double x = xy[0], y = xy[1];
				int errcode = transform_point(from,to,x,y);
				out[0] = x; out[1] = y;
				return errcode;
			}

			inline int transform_point(const char *from, const std::string &to, double xy[], double out[]) { 
				double x = xy[0], y = xy[1];
				int errcode = transform_point(from,to,x,y);
				out[0] = x; out[1] = y;
				return errcode;
			}

			inline int transform_point(const std::string &from, const char *to, double xy[], double out[]) { 
				double x = xy[0], y = xy[1];
				int errcode = transform_point(from,to,x,y);
				out[0] = x; out[1] = y;
				return errcode;
			}

			inline int transform_point(const std::string &from, const std::string &to, double xy[], double out[]) { 
				double x = xy[0], y = xy[1];
				int errcode = transform_point(from,to,x,y);
				out[0] = x; out[1] = y;
				return errcode;
			}

			inline int transform_point(const char *from, const char *to, double &x, double &y) { 
				return proj_transform_point(projs[std::string(from)],projs[std::string(to)],&x,&y);
			}

			inline int transform_point(const char *from, const std::string &to, double &x, double &y) { 
				return proj_transform_point(projs[std::string(from)],projs[to],&x,&y);
			}

			inline int transform_point(const std::string &from, const char *to, double &x, double &y) { 
				return proj_transform_point(projs[from],projs[std::string(to)],&x,&y);
			}

			inline int transform_point(const std::string &from, const std::string &to, double &x, double &y) { 
				return proj_transform_point(projs[from],projs[to],&x,&y);
			}

			inline int proj_transform_point(const char *from, const char *to, double *x, double *y) {
				if (from==NULL || to==NULL) return PROJ_ERR;
				
				source = pj_init_plus(from); 
				target = pj_init_plus(to);
				alloc = true;
				if (source==NULL || target==NULL) return PROJ_ERR;

				int retval = pj_transform(source, target, 1, 1, x, y, NULL); // 0 on success
				free();
				return (retval == 0) ? PROJ_OK : PROJ_ERR;
			}

		private:
			bool alloc;
			projPJ source, target;

			// map containing all the default projections from the Suite
			std::map<std::string,const char*> projs = {
				{"degrees"     ,"+ellps=WGS84 +a=57.29577951308232 +proj=eqc +lon_0=0.0 +no_defs"},                                          // Plate Carree
				{"mercator"    ,"+ellps=WGS84 +proj=merc +lon_0=0.0 +x_0=989634.3811336625 +y_0=-3512473.95569 +units=m +no_defs"},          // Mercator centered on MED
				{"cylindrical" ,"+datum=WGS84 +ellps=WGS84 +proj=eqc +lat_ts=0 +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +units=m +no_defs +no_defs"}, // Cylindrical using EPSG:4087
				{"google"      ,"+a=6378137.0 +b=6378137.0 +nadgrids=@null +proj=merc +lon_0=0.0 +x_0=0.0 +y_0=0.0 +units=m +no_defs"},      // Google Mercator
				{"mollweide"   ,"+a=6378137.0 +proj=moll +lon_0=0 +no_defs"},
				{"orthographic","+ellps=WGS84 +proj=ortho +lon_0=0.0 +lat_0=0.0 +no_defs"},
				{"robinson"    ,"+a=6378137.0 +proj=robin +lon_0=0 +no_defs"},
				{"satellite"   ,"+a=6378137.0 +proj=nsper +lon_0=17.5 +lat_0=36.4 +h=6779000 +x_0=0 +y_0=0 +units=m +no_defs"},              // Nearside Perspective centered on MED at ISS altitude 
				{"eckert iv"   ,"+a=6378137.0 +proj=eck4 +lon_0=0 +no_defs"},
				{"equal earth" ,"+ellps=WGS84 +proj=eqearth +lon_0=0 +no_defs"},
				{"epsg 3857"   ,"+ellps=WGS84 +a=6378137 +b=6378137 +nadgrids=@null +proj=merc +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +wktext + +no_defs +no_defs"}
			};

			inline void free() { if (alloc) {pj_free(source); pj_free(target); alloc=false;} }
	};
}

#endif