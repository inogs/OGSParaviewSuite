/*=========================================================================

  Module:    OGS Main class

  Library to deal with the reading and writing of OGS files as well as the
  NETCDF4 interface. 

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef OGS_H
#define OGS_H

#include <vector>
#include <string>
#include "field.h"

/* SIZES FOR MALLOCS */

#define CHRSZ sizeof(char)
#define INTSZ sizeof(int)
#define DBLSZ sizeof(double)

/* BUFFER SIZE FOR FILE READING */

#define VARSZ   80
#define BUFFSZ  256
#define LBUFFSZ 512

/* MACRO FOR ERROR */

//#define ERROR(errstr) {fprintf(stderr,"Error in %s line %d:\n%s\n",__FILE__,__LINE__,errstr); exit(-1);}
#define ERROR(errstr,errid) {std::fprintf(stderr,"Error in %s line %d:\n%s\n",__FILE__,__LINE__,errstr); return(-errid);}

namespace ogs
{
	/* OGS VARIABLE

		Stores an array with the name of the variable, its name inside
		the NETCDF and the relative path to the netcdf file.

	*/
	class OGS_VAR {
		public:
			inline OGS_VAR()                         {};
			inline ~OGS_VAR()                        {};

			inline int  get_nvars()                  { return this->_n; }
			inline void set_nvars(int n)             { this->_n = n; this->allocate(); }

			inline void  set_name(int i, char* str)  { this->_name[i]  = std::string(str); }
			inline void  set_vname(int i, char* str) { this->_vname[i] = std::string(str); }
			inline void  set_path(int i, char* str)  { this->_path[i]  = std::string(str); }

			inline const char *get_name(int i)       { return this->_name[i].c_str();  }
			inline const char *get_vname(int i)      { return this->_vname[i].c_str(); }
			inline const char *get_path(int i)       { return this->_path[i].c_str();  }

		private:
			int _n = 0;
			std::vector<std::string> _name, _vname, _path;

			inline void allocate() { this->_name.resize(this->_n); this->_vname.resize(this->_n); this->_path.resize(this->_n); }
	};

	/* OGS CLASS

		Stores the necessary information from the OGS master file and helps in reading/writing 
		the variables.

	*/
	class OGS {
		public:
			// Constructors
			inline OGS();
			inline OGS(const char *fname);

			// Destructor
			inline ~OGS();

			// Set/Get methods
			inline void SetFile(const char *fname);
			inline void SetWdir(const char *fname);
			inline void SetMfile(const char *fname);
			inline void Setlon2m(const int n, double *arr);
			inline void Setlat2m(const int n, double *arr);
			inline void Setnavlev(const int n, double *arr);
			inline void SetMask(const int i, double *mask);
			inline void Setncells();

			inline std::string meshfile();
			inline std::string meshmask();

			inline int  nlon();
			inline int  nlat();
			inline int  nlev();
			inline int  ncells();

			inline double  lon2meters(int i);
			inline double *lon2meters();
			inline double  lat2meters(int i);
			inline double *lat2meters();
			inline double  nav_lev(int i);
			inline double *nav_lev();

			inline field::Field<double> &mask(int i);
			inline field::Field<double> &e1();
			inline field::Field<double> &e2();
			inline field::Field<double> &e3();

			inline int         var_n(int i);
			inline const char *var_name(int i, int j);
			inline const char *var_vname(int i, int j);
			inline std::string var_path(int i, int j, int t);

			inline int         ntsteps();
			inline const char *datetime(int i);

			// Functions
			int  readMainFile();
			int  readMesh();
			int  writeMesh();
			void readMeshmask();

			void print();

		private:
			std::string _ogsfile, _wrkdir, _meshfile, _meshmask;
			OGS_VAR _vars[4]; // Variables information
			int _ntsteps;
			std::vector<std::string> _datetime;

			// Mesh size information
			int _nlat, _nlon, _nlev, _ncells;
			std::vector<double> _lon2m, _lat2m, _nav_lev;
			field::Field<double> _masks[2];
			field::Field<double> _e1, _e2, _e3;

			std::string var_WritePath(int i, int j, const char *str, const char *token);
	};

	// Constructions and destructors
	inline         OGS::OGS()                           {}
	inline         OGS::OGS(const char *fname)          { this->SetFile(fname); }
	inline         OGS::~OGS()                          {}

	// Set/Get methods
	inline void    OGS::SetFile(const char *fname) { this->_ogsfile = std::string(fname); }
	inline void    OGS::SetWdir(const char *fname) { this->_wrkdir = std::string(fname); }
	inline void    OGS::SetMfile(const char *fname){ this->_meshfile = std::string(fname); }
	inline void    OGS::Setncells()                { this->_ncells = (this->_nlon-1)*(this->_nlat-1)*(this->_nlev-1); }
	inline std::string OGS::meshfile()             { return (this->_wrkdir + std::string("/") + this->_meshfile); }
	inline std::string OGS::meshmask()             { return (this->_wrkdir + std::string("/") + this->_meshmask); }
	inline int     OGS::nlon()                     { return this->_nlon; }
	inline int     OGS::nlat()                     { return this->_nlat; }
	inline int     OGS::nlev()                     { return this->_nlev; }
	inline int     OGS::ncells()                   { return this->_ncells; }
	inline double  OGS::lon2meters(int i)          { return (i >= 0) ? this->_lon2m[i] : this->_lon2m[this->_nlon + i]; }
	inline double *OGS::lon2meters()               { return this->_lon2m.data(); }
	inline double  OGS::lat2meters(int i)          { return (i >= 0) ? this->_lat2m[i] : this->_lat2m[this->_nlat + i];}
	inline double *OGS::lat2meters()               { return this->_lat2m.data(); }
	inline double  OGS::nav_lev(int i)             { return (i >= 0) ? this->_nav_lev[i] : this->_nav_lev[this->_nlev + i];}
	inline double *OGS::nav_lev()                  { return this->_nav_lev.data(); }
	inline field::Field<double> &OGS::mask(int i)  { return this->_masks[i]; }
	inline field::Field<double> &OGS::e1()         { return this->_e1; }
	inline field::Field<double> &OGS::e2()         { return this->_e2; }
	inline field::Field<double> &OGS::e3()         { return this->_e3; }
	inline int     OGS::var_n(int i)               { return this->_vars[i].get_nvars(); }
	inline const char *OGS::var_name(int i, int j) { return this->_vars[i].get_name(j); }
	inline const char *OGS::var_vname(int i, int j){ return this->_vars[i].get_vname(j); }
	inline int     OGS::ntsteps()                  { return this->_ntsteps; }
	inline const char *OGS::datetime(int i)        { return this->_datetime[i].c_str(); }
	inline std::string OGS::var_path(int i, int j, int t) { 
		return(this->_wrkdir + std::string("/") + this->var_WritePath(i,j,this->datetime(t),"*"));
	}
	inline void    OGS::Setlon2m(const int n, double *arr) {
		this->_nlon = n;
		this->_lon2m.resize(n); 
		if (arr != NULL) 
			std::memcpy(this->_lon2m.data(),arr,n*sizeof(double));
	}
	inline void    OGS::Setlat2m(const int n, double *arr) {
		this->_nlat = n;
		this->_lat2m.resize(n); 
		if (arr != NULL) 
			std::memcpy(this->_lat2m.data(),arr,n*sizeof(double));
	}
	inline void    OGS::Setnavlev(const int n, double *arr) {
		this->_nlev = n;
		this->_nav_lev.resize(n); 
		if (arr != NULL) 
			std::memcpy(this->_nav_lev.data(),arr,n*sizeof(double));
	}
	inline void    OGS::SetMask(const int i, double *mask) {
		this->_masks[i].set_dim(this->_ncells,1);
		std::memcpy(this->_masks[i].data(),mask,this->_ncells*sizeof(double));
	}
}
#endif