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

#include <cstdint>
#include <vector>
#include <algorithm>
#include <string>

#include "field.h"

/* SIZES FOR MALLOCS */

#define CHRSZ sizeof(char)
#define INTSZ sizeof(int)
#define UI8SZ sizeof(uint8_t)
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

			inline void  set_name(int i, char *str)  { this->_name[i]  = std::string(str); }
			inline void  set_vname(int i, char *str) { this->_vname[i] = std::string(str); }
			inline void  set_path(int i, char *str)  { this->_path[i]  = std::string(str); }

			inline const char *get_name(int i)       { return this->_name[i].c_str();  }
			inline const char *get_vname(int i)      { return this->_vname[i].c_str(); }
			inline const char *get_path(int i)       { return this->_path[i].c_str();  }

			inline int find_name(const char *str) {
				for (int ii = 0; ii < this->_name.size(); ++ii)
					if (this->_name[ii] == std::string(str)) 
						return  ii;
				return -1;
			}

		private:
			int _n = 0;
			std::vector<std::string> _name, _vname, _path;

			inline void allocate() { this->_name.resize(this->_n); this->_vname.resize(this->_n); this->_path.resize(this->_n); }
	};

	/* OGS MESH

		Stores information on the mesh and the projection used

	*/
	class OGS_MESH {
		public:
			inline OGS_MESH()                         {};
			inline ~OGS_MESH()                        {};

			inline void set_name(const char *str)     { this->_name     = std::string(str); }
			inline void set_meshfile(const char *str) { this->_meshfile = std::string(str); }
			inline void set_meshmask(const char *str) { this->_meshmask = std::string(str); }

			inline std::string get_name()             { return this->_name;     }
			inline std::string get_meshfile()         { return this->_meshfile; }
			inline std::string get_meshmask()         { return this->_meshmask; }

		private:
			std::string _name, _meshfile, _meshmask;
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
			inline OGS(std::string fname);

			// Destructor
			inline ~OGS();

			// Set/Get methods
			inline void SetFile(const char *fname);
			inline void SetWdir(const char *fname);
			inline void SetMfile(const char *fname);
			inline void Setlon2m(const int n, double *arr);
			inline void Setlat2m(const int n, double *arr);
			inline void Setnavlev(const int n, double *arr);
			inline void SetMask(const int i, const int m, uint8_t *mask);
			inline void Setncells();

			inline std::string projection(const int i);
			inline std::string meshfile(const int i);
			inline std::string meshmask(const int i);

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

			inline field::Field<uint8_t> &mask(int i);
			inline field::Field<double> &e1();
			inline field::Field<double> &e2();
			inline field::Field<double> &e3();

			inline int         var_n(int i);
			inline const char *var_name(int i, int j);
			inline const char *var_vname(int i, int j);
			inline const char *var_vname(const char *name);
			inline const char *var_vname(std::string name);
			inline std::string var_path(int i, int j, int t);
			inline std::string var_path(const char *vname, int t);
			inline std::string var_path(std::string vname, int t);

			inline int         ntsteps();
			inline const char *datetime(int i);

			// Functions
			int  readMainFile();
			int  readMesh(const int i);
			int  writeMesh(const int i);
			void readMeshmask(const int i);

			void print();

		private:
			std::string _ogsfile, _wrkdir;
			OGS_VAR _vars[4]; // Variables information
			int _ntsteps;
			std::vector<std::string> _datetime;
			std::vector<OGS_MESH> mesh_data;

			// Mesh size information
			int _nlat, _nlon, _nlev, _ncells;
			std::vector<double> _lon2m, _lat2m, _nav_lev;
			field::Field<uint8_t> _masks[2];
			field::Field<double> _e1, _e2, _e3;

			std::string var_WritePath(int i, int j, const char *str, const char *token);
			std::string var_WritePath(const char *vname, const char *str, const char *token);
	};

	// Constructions and destructors
	inline         OGS::OGS()                           {}
	inline         OGS::OGS(const char *fname)          { this->SetFile(fname); }
	inline         OGS::OGS(std::string fname)          { this->SetFile(fname.c_str()); }
	inline         OGS::~OGS()                          {}

	// Set/Get methods
	inline void    OGS::SetFile(const char *fname) { this->_ogsfile = std::string(fname); }
	inline void    OGS::SetWdir(const char *fname) { this->_wrkdir = std::string(fname); }
	inline void    OGS::SetMfile(const char *fname){ this->mesh_data.resize(1); this->mesh_data[0].set_meshfile(fname); }
	inline void    OGS::Setncells()                { this->_ncells = (this->_nlon-1)*(this->_nlat-1)*(this->_nlev-1); }
	inline std::string OGS::projection(const int i){ return mesh_data[i].get_name(); }
	inline std::string OGS::meshfile(const int i)  { return (this->_wrkdir + std::string("/") + mesh_data[i].get_meshfile()); }
	inline std::string OGS::meshmask(const int i)  { return (this->_wrkdir + std::string("/") + mesh_data[i].get_meshmask()); }
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
	inline field::Field<uint8_t> &OGS::mask(int i) { return this->_masks[i]; }
	inline field::Field<double>  &OGS::e1()        { return this->_e1; }
	inline field::Field<double>  &OGS::e2()        { return this->_e2; }
	inline field::Field<double>  &OGS::e3()        { return this->_e3; }
	inline int         OGS::ntsteps()              { return this->_ntsteps; }
	inline const char *OGS::datetime(int i)        { return this->_datetime[i].c_str(); }
	inline int         OGS::var_n(int i)           { return this->_vars[i].get_nvars(); }
	inline const char *OGS::var_name(int i, int j) { return this->_vars[i].get_name(j); }
	inline const char *OGS::var_vname(int i, int j){ return this->_vars[i].get_vname(j); }
	inline const char *OGS::var_vname(const char *name) { 
		for (int i = 0; i < 4; ++i) {
			int j = this->_vars[i].find_name(name);
			if (j >= 0) return this->_vars[i].get_vname(j);
		}
		return NULL; 
	}
	inline const char *OGS::var_vname(std::string name) { return this->var_vname(name.c_str()); }
	inline std::string OGS::var_path(int i, int j, int t) { 
		return(this->_wrkdir + std::string("/") + this->var_WritePath(i,j,this->datetime(t),"*"));
	}
	inline std::string OGS::var_path(const char *name, int t) { 
		return(this->_wrkdir + std::string("/") + this->var_WritePath(name,this->datetime(t),"*"));
	}
	inline std::string OGS::var_path(std::string name, int t) { return this->var_path(name.c_str(),t); }
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
	inline void    OGS::SetMask(const int i, const int m, uint8_t *mask) { this->_masks[i].set_dim(this->_ncells,m); this->_masks[i].set_val(mask); }
}
#endif