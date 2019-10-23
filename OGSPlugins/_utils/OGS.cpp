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

#include <cstdlib>
#include <cstdio>
#include <cstring>

#include <vector>
#include <string>

#ifdef __linux__
// Include OpenMP when working with GCC
#include <omp.h>
#endif

#include "OGS.hpp"

#ifndef OGS_NO_NETCDF
#include "netcdfio.hpp"
#endif

#define CLLIND(ii,jj,kk,nx,ny) ( (nx-1)*(ny-1)*(kk) + (nx-1)*(jj) + (ii) )
#define PNTIND(ii,jj,kk,nx,ny) ( (nx)*(ny)*(kk) + (nx)*(jj) + (ii) )

using namespace ogs;

/* PROTOTYPES */

char *reads(char *line, int size, FILE *fin);
char *trim(char *str);

/* CLASS FUNCTIONS */

int OGS::readMainFile() {

	char  line[BUFFSZ];
	char *linep;

	int sec = 0, n_var_r = -1, n_time_r = -1;

	bool bworkdir = false, bmesh = false, btime = false;
	bool bavephys = false, bavefreq = false, bforcings = false, bgenerals = false;

	// Open file for reading
	FILE *myfile;
	myfile = std::fopen(this->_ogsfile.c_str(),"r"); if (myfile == NULL) ERROR("Cannot open file",1)

	// Read file line by line
	while(reads(line,sizeof(line),myfile)) {
		
		/* IGNORE COMMENTS */
		if (line[0] == '#') continue;

		/* DEFINITION OF SECTIONS */
		if (std::string(line) == "WRKDIR")   {sec++; continue;} // Section containing the working directory
		if (std::string(line) == "MESH")     {sec++; continue;} // Section containing the mesh information
		if (std::string(line) == "AVE_PHYS") {sec++; continue;} // Section containing the physical variables
		if (std::string(line) == "AVE_FREQ") {sec++; continue;} // Section containing the biogeochemical variables
		if (std::string(line) == "FORCINGS") {sec++; continue;} // Section containing the forcings variables
		if (std::string(line) == "GENERALS") {sec++; continue;} // Section containing general variables (e.g., downloaded from web)
		if (std::string(line) == "TIME")     {sec++; continue;} // Section containing the time-stepping information

		/* WORKDIR SECTION */
		if (sec == 1 && !bworkdir) { this->_wrkdir = std::string(trim(line)); bworkdir = true; }

		/* MESH SECTION */
		if (sec == 2 && !bmesh) {
			if (n_var_r >= 0) {
				if (n_var_r > this->mesh_data.size() - 1) { bmesh = true; n_var_r = -1; continue; }
				// Read projection name
				strtok(line,":");         this->mesh_data[n_var_r].set_name(trim(line));
				// Read mesh file
				linep = strtok(NULL,":"); this->mesh_data[n_var_r].set_meshfile(trim(linep));
				// Read meshmask file
				linep = strtok(NULL," "); this->mesh_data[n_var_r].set_meshmask(trim(linep));
				n_var_r++;
			} else {
				this->mesh_data.resize( std::atoi(line) );
				n_var_r++;
			}
		}

		/* AVE PHYS SECTION */
		if (sec == 3 && !bavephys) {
			if (n_var_r >= 0) {
				if (n_var_r > this->_vars[0].get_nvars() - 1) { bavephys = true; n_var_r = -1; continue; }
				// Split the string, read variable name
				strtok(line,":");         this->_vars[0].set_name(n_var_r,trim(line));
				// Split the string, read netcfd name
				linep = strtok(NULL,":"); this->_vars[0].set_vname(n_var_r,trim(linep));
				// Split the string, read variable path
				linep = strtok(NULL," "); this->_vars[0].set_path(n_var_r,trim(linep));
				n_var_r++;
			} else {
				this->_vars[0].set_nvars( std::atoi(line) );
				n_var_r++;
			}
		}

		/* AVE FREQ SECTION */
		if (sec == 4 && !bavefreq) {
			if (n_var_r >= 0) {
				if (n_var_r > this->_vars[1].get_nvars() - 1) { bavefreq = true; n_var_r = -1; continue; }
				// Split the string, read variable name
				strtok(line,":");         this->_vars[1].set_name(n_var_r,trim(line));
				// Split the string, read netcfd name
				linep = strtok(NULL,":"); this->_vars[1].set_vname(n_var_r,trim(linep));
				// Split the string, read variable path
				linep = strtok(NULL," "); this->_vars[1].set_path(n_var_r,trim(linep));
				n_var_r++;
			} else {
				this->_vars[1].set_nvars( std::atoi(line) );
				n_var_r++;
			}
		}

		/* FORCINGS SECTION */
		if (sec == 5 && !bforcings) {
			if (n_var_r >= 0) {
				if (n_var_r > this->_vars[2].get_nvars() - 1) { bforcings = true; n_var_r = -1; continue; }
				// Split the string, read variable name
				strtok(line,":");         this->_vars[2].set_name(n_var_r,trim(line));
				// Split the string, read netcfd name
				linep = strtok(NULL,":"); this->_vars[2].set_vname(n_var_r,trim(linep));
				// Split the string, read variable path
				linep = strtok(NULL," "); this->_vars[2].set_path(n_var_r,trim(linep));
				n_var_r++;
			} else {
				this->_vars[2].set_nvars( std::atoi(line) );
				n_var_r++;
			}
		}

		/* GENERALS SECTION */
		if (sec == 6 && !bgenerals) {
			if (n_var_r >= 0) {
				if (n_var_r > this->_vars[3].get_nvars() - 1) { bgenerals = true; n_var_r = -1; continue; }
				// Split the string, read variable name
				strtok(line,":");         this->_vars[3].set_name(n_var_r,trim(line));
				// Split the string, read netcfd name
				linep = strtok(NULL,":"); this->_vars[3].set_vname(n_var_r,trim(linep));
				// Split the string, read variable path
				linep = strtok(NULL," "); this->_vars[3].set_path(n_var_r,trim(linep));
				n_var_r++;
			} else {
				this->_vars[3].set_nvars( std::atoi(line) );
				n_var_r++;
			}
		}

		/* TIMESTEP SECTION */
		if (sec == 7 && !btime) {
			if (n_time_r >= 0) {
				if (n_time_r > this->_ntsteps - 1) { btime = true; continue; }
				// Read datetime
				this->_datetime[n_time_r] = std::string(trim(line));
				n_time_r++;
			} else {
				// Read number of timesteps
				this->_ntsteps = std::atoi(line);
				// Allocate
				this->_datetime.resize(this->_ntsteps);
				n_time_r++;
			}
		}
	}

	std::fclose(myfile);
	return 1;
}

int OGS::readMesh(const int i) {
	// Open file for reading
	FILE *myfile;
	myfile = std::fopen(this->meshfile(i).c_str(),"rb"); if (myfile == NULL) ERROR("Cannot open file.",1)

	// Read the dimensions
	if (std::fread(&this->_nlon,INTSZ,1,myfile) != 1) ERROR("Error reading file",2) 
	if (std::fread(&this->_nlat,INTSZ,1,myfile) != 1) ERROR("Error reading file",2) 
	if (std::fread(&this->_nlev,INTSZ,1,myfile) != 1) ERROR("Error reading file",2)
	this->Setncells();

	// Read the vectors
	this->_lon2m.resize(this->_nlon,0);
	if (std::fread(this->_lon2m.data(),DBLSZ,this->_nlon,myfile)   != this->_nlon) ERROR("Error reading file",2)
	this->_lat2m.resize(this->_nlat,0);
	if (std::fread(this->_lat2m.data(),DBLSZ,this->_nlat,myfile)   != this->_nlat) ERROR("Error reading file",2)
	this->_nav_lev.resize(this->_nlev,0);	
	if (std::fread(this->_nav_lev.data(),DBLSZ,this->_nlev,myfile) != this->_nlev) ERROR("Error reading file",2)

	// Read the masks
	for (int ii = 0; ii < 3; ++ii) {
		int m;
		if (std::fread(&m,INTSZ,1,myfile) != 1) ERROR("Error reading file",2)
		this->_masks[ii].set_dim(this->_ncells,m);
		if (std::fread(this->_masks[ii].data(),UI8SZ,this->_masks[ii].get_sz(),myfile) != this->_masks[ii].get_sz()) ERROR("Error reading file",2)
	}

	std::fclose(myfile);
	return 1;
}

int OGS::writeMesh(const int i) {
	// Open file for writing
	FILE *myfile;
	myfile = std::fopen(this->meshfile(i).c_str(),"wb"); if (myfile == NULL) ERROR("Cannot open file.",1)

	// Write the dimensions
	if (std::fwrite(&this->_nlon,INTSZ,1,myfile) != 1) ERROR("Error writing file",2) 
	if (std::fwrite(&this->_nlat,INTSZ,1,myfile) != 1) ERROR("Error writing file",2) 
	if (std::fwrite(&this->_nlev,INTSZ,1,myfile) != 1) ERROR("Error writing file",2)

	// Write the vectors
	if (std::fwrite(this->_lon2m.data(),DBLSZ,this->_nlon,myfile)   != this->_nlon) ERROR("Error writing file",2)
	if (std::fwrite(this->_lat2m.data(),DBLSZ,this->_nlat,myfile)   != this->_nlat) ERROR("Error writing file",2)
	if (std::fwrite(this->_nav_lev.data(),DBLSZ,this->_nlev,myfile) != this->_nlev) ERROR("Error writing file",2)

	// Write the masks
	for (int ii = 0; ii < 3; ++ii) {
		int m = this->_masks[ii].get_m();
		if (std::fwrite(&m,INTSZ,1,myfile) != 1) ERROR("Error writing file",2)
		if (std::fwrite(this->_masks[ii].data(),UI8SZ,this->_masks[ii].get_sz(),myfile) != this->_masks[ii].get_sz()) ERROR("Error writing file",2)
	}

	std::fclose(myfile);
	return 1;
}

void OGS::readMeshmask(const int i) {
#ifndef OGS_NO_NETCDF
	// Read the fields as double arrays
	int ncells2D = 1*(this->_nlon-1)*(this->_nlat-1);
	
	std::vector<double> e1t(ncells2D), e1u(ncells2D), e1v(ncells2D), e1f(ncells2D);
	NetCDF::readNetCDF(this->meshmask(i).c_str(),"e1t",e1t);
	NetCDF::readNetCDF(this->meshmask(i).c_str(),"e1u",e1u);
	NetCDF::readNetCDF(this->meshmask(i).c_str(),"e1v",e1v);
	NetCDF::readNetCDF(this->meshmask(i).c_str(),"e1f",e1f);

	std::vector<double> e2t(ncells2D), e2u(ncells2D), e2v(ncells2D), e2f(ncells2D);
	NetCDF::readNetCDF(this->meshmask(i).c_str(),"e2t",e2t);
	NetCDF::readNetCDF(this->meshmask(i).c_str(),"e2u",e2u);
	NetCDF::readNetCDF(this->meshmask(i).c_str(),"e2v",e2v);
	NetCDF::readNetCDF(this->meshmask(i).c_str(),"e2f",e2f);

	std::vector<double> e3t(this->_ncells), e3u(this->_ncells), e3v(this->_ncells), e3w(this->_ncells);
	NetCDF::readNetCDF(this->meshmask(i).c_str(),"e3t_0",e3t);
	NetCDF::readNetCDF(this->meshmask(i).c_str(),"e3u_0",e3u);
	NetCDF::readNetCDF(this->meshmask(i).c_str(),"e3v_0",e3v);
	NetCDF::readNetCDF(this->meshmask(i).c_str(),"e3w_0",e3w);

	// Allocate the fields
	this->_e1.set_dim(this->_ncells,4);
	this->_e2.set_dim(this->_ncells,4);
	this->_e3.set_dim(this->_ncells,4);

	// Fill the fields
	#pragma omp parallel for collapse(3)
	for (int kk = 0; kk < this->_nlev-1; ++kk) {
		for (int jj = 0; jj < this->_nlat-1; ++jj) {
			for (int ii = 0; ii < this->_nlon-1; ++ii) {
				// e1
				this->_e1[CLLIND(ii,jj,kk,this->_nlon,this->_nlat)][0] = e1t[CLLIND(ii,jj,0,this->_nlon,this->_nlat)];
				this->_e1[CLLIND(ii,jj,kk,this->_nlon,this->_nlat)][1] = e1u[CLLIND(ii,jj,0,this->_nlon,this->_nlat)];
				this->_e1[CLLIND(ii,jj,kk,this->_nlon,this->_nlat)][2] = e1v[CLLIND(ii,jj,0,this->_nlon,this->_nlat)];
				this->_e1[CLLIND(ii,jj,kk,this->_nlon,this->_nlat)][3] = e1f[CLLIND(ii,jj,0,this->_nlon,this->_nlat)];
				// e2
				this->_e2[CLLIND(ii,jj,kk,this->_nlon,this->_nlat)][0] = e2t[CLLIND(ii,jj,0,this->_nlon,this->_nlat)];
				this->_e2[CLLIND(ii,jj,kk,this->_nlon,this->_nlat)][1] = e2u[CLLIND(ii,jj,0,this->_nlon,this->_nlat)];
				this->_e2[CLLIND(ii,jj,kk,this->_nlon,this->_nlat)][2] = e2v[CLLIND(ii,jj,0,this->_nlon,this->_nlat)];
				this->_e2[CLLIND(ii,jj,kk,this->_nlon,this->_nlat)][3] = e2f[CLLIND(ii,jj,0,this->_nlon,this->_nlat)];
				// e3
				this->_e3[CLLIND(ii,jj,kk,this->_nlon,this->_nlat)][0] = e3t[CLLIND(ii,jj,kk,this->_nlon,this->_nlat)];
				this->_e3[CLLIND(ii,jj,kk,this->_nlon,this->_nlat)][1] = e3u[CLLIND(ii,jj,kk,this->_nlon,this->_nlat)];
				this->_e3[CLLIND(ii,jj,kk,this->_nlon,this->_nlat)][2] = e3v[CLLIND(ii,jj,kk,this->_nlon,this->_nlat)];
				this->_e3[CLLIND(ii,jj,kk,this->_nlon,this->_nlat)][3] = e3w[CLLIND(ii,jj,kk,this->_nlon,this->_nlat)];
			}
		}
	}
#endif
}

std::string OGS::var_WritePath(int i, int j, const char *str, const char *token) {
	std::string aux = std::string(this->_vars[i].get_path(j));
	int pos = aux.find(token);

	return aux.replace(pos,1,str);
}

std::string OGS::var_WritePath(const char *vname, const char *str, const char *token) {
	int i, j;

	for (i = 0; i < 4; ++i) {
		j = this->_vars[i].find_name(vname);
		if (j >= 0) break;
	}

	std::string aux = std::string(this->_vars[i].get_path(j));
	int pos = aux.find(token);

	return aux.replace(pos,1,str);
}

void OGS::print() {
	std::printf("OGS CLASS\n");
	std::printf("---------\n");
	std::printf("Main file:          %s\n", this->_ogsfile.c_str());
	std::printf("Working directory:  %s\n", this->_wrkdir.c_str());

	std::printf("Mesh information\n");
	for (OGS_MESH mesh : this->mesh_data) {
		std::printf("    Projection:     %s\n", mesh.get_name().c_str());
		std::printf("    Mesh file:      %s\n", mesh.get_meshfile().c_str());
		std::printf("    Meshmask file:  %s\n", mesh.get_meshmask().c_str());
	}
	this->readMesh(0);
	std::printf("    Dimensions:     (%d,%d,%d)\n", this->_nlon, this->_nlat, this->_nlev);
	std::printf("    Longitude [m]:  (%.2f,%.2f)\n", this->lon2meters(0), this->lon2meters(-1));
	std::printf("    Latitude [m]:   (%.2f,%.2f)\n", this->lat2meters(0), this->lat2meters(-1));
	std::printf("    Depth [m]:      (%.2f,%.2f)\n", this->nav_lev(0),    this->nav_lev(-1));

	std::printf("Variables information\n");
	for (int ii = 0; ii < 4; ii++) {
		std::printf("    Variables %d:    %d\n",ii,this->var_n(ii));
		for (int jj = 0; jj < this->var_n(ii); jj++)
			std::printf("      %d Name:       %s (%s)\n",jj,this->var_name(ii,jj),this->var_vname(ii,jj));
	}

}

/* INTERFACE C FUNCTIONS FOR PYTHON LIBRARY */

extern "C"
{
	int OGSWriteMesh(const char *fname, const char *wrkdir, const int nlon, const int nlat, const int nlev,
		double *lon2m, double *lat2m, double *nav_lev, uint8_t *bmask, uint8_t *cmask, uint8_t *lmask) {
		// Create a new instance of the class
		OGS ogscls;
		// Populate the class
		ogscls.SetMfile(fname);
		ogscls.SetWdir(wrkdir);
		ogscls.Setlon2m(nlon,lon2m);
		ogscls.Setlat2m(nlat,lat2m);
		ogscls.Setnavlev(nlev,nav_lev);
		ogscls.Setncells();
		// Load the masks
		ogscls.SetMask(0,16,bmask);
		ogscls.SetMask(1,1,cmask);
		ogscls.SetMask(2,1,lmask);
		// Print
		return ogscls.writeMesh(0);
	}
}

/* AUXILIARY FUNCTIONS */

char *reads(char *line, int size, FILE *fin) {
	char *l = std::fgets(line,size,fin);
	if (l != NULL)
	  for (int ii = 0; line[ii] != '\0'; ii++)
	    if (line[ii] == 10) { line[ii] = '\0'; break; }
	return l;
}
char *trim(char *str) {
        char *end;
        while(std::isspace(*str)) str++; // Trim leading space
        if(*str == 0) return str;// All spaces?
        // Trim trailing space
        end = str + std::strlen(str) - 1;
        while(end > str && std::isspace(*end)) end--;
        // Write new null terminator
        *(end+1) = 0;
        return str;
}