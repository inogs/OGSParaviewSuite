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

#include <string>

#include "netcdfio.hpp"
#include "OGS.hpp"

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
	myfile = std::fopen(this->_ogsfile.c_str(),"r"); if (myfile == NULL) ERROR2("Cannot open file",1)

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
			// Read mesh file
			strtok(line,":");         this->_meshfile = std::string(trim(line));
			// Read meshmask file
			linep = strtok(NULL," "); this->_meshmask = std::string(trim(linep));
			bmesh = true;
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

int OGS::readMesh() {
	// Open file for reading
	FILE *myfile;
	myfile = std::fopen(this->meshfile().c_str(),"rb"); if (myfile == NULL) ERROR2("Cannot open file.",1)

	// Read the dimensions
	if (std::fread(&this->_nlon,INTSZ,1,myfile) != 1) ERROR2("Error reading file",2) 
	if (std::fread(&this->_nlat,INTSZ,1,myfile) != 1) ERROR2("Error reading file",2) 
	if (std::fread(&this->_nlev,INTSZ,1,myfile) != 1) ERROR2("Error reading file",2)
	this->Setncells();

	// Read the vectors
	this->_lon2m.resize(this->_nlon,0);
	if (std::fread(this->_lon2m.data(),DBLSZ,this->_nlon,myfile)   != this->_nlon) ERROR2("Error reading file",2)
	this->_lat2m.resize(this->_nlat,0);
	if (std::fread(this->_lat2m.data(),DBLSZ,this->_nlat,myfile)   != this->_nlat) ERROR2("Error reading file",2)
	this->_nav_lev.resize(this->_nlev,0);	
	if (std::fread(this->_nav_lev.data(),DBLSZ,this->_nlev,myfile) != this->_nlev) ERROR2("Error reading file",2)

	// Read the masks
	for (int ii = 0; ii < 2; ii++) {
		this->_masks[ii].set_dim(this->_ncells,1);
		if (std::fread(this->_masks[ii].data(),DBLSZ,this->_ncells,myfile) != this->_ncells) ERROR2("Error reading file",2)
	}

	std::fclose(myfile);
	return 1;
}

int OGS::writeMesh() {
	// Open file for writing
	FILE *myfile;
	myfile = std::fopen(this->meshfile().c_str(),"wb"); if (myfile == NULL) ERROR2("Cannot open file.",1)

	// Write the dimensions
	if (std::fwrite(&this->_nlon,INTSZ,1,myfile) != 1) ERROR2("Error writing file",2) 
	if (std::fwrite(&this->_nlat,INTSZ,1,myfile) != 1) ERROR2("Error writing file",2) 
	if (std::fwrite(&this->_nlev,INTSZ,1,myfile) != 1) ERROR2("Error writing file",2)

	// Write the vectors
	if (std::fwrite(this->_lon2m.data(),DBLSZ,this->_nlon,myfile)   != this->_nlon) ERROR2("Error writing file",2)
	if (std::fwrite(this->_lat2m.data(),DBLSZ,this->_nlat,myfile)   != this->_nlat) ERROR2("Error writing file",2)
	if (std::fwrite(this->_nav_lev.data(),DBLSZ,this->_nlev,myfile) != this->_nlev) ERROR2("Error writing file",2)

	// Write the masks
	for (int ii = 0; ii < 2; ii++)
		if (std::fwrite(this->_masks[ii].data(),DBLSZ,this->_ncells,myfile) != this->_ncells) ERROR2("Error writing file",2)

	std::fclose(myfile);
	return 1;
}

void OGS::readMeshmask() {
	// Read the fields as double arrays
	double *e1t, *e1u, *e1v, *e1f;
	e1t = NetCDF::readNetCDF(this->meshmask().c_str(),"e1t",1*(this->_nlon-1)*(this->_nlat-1));
	e1u = NetCDF::readNetCDF(this->meshmask().c_str(),"e1u",1*(this->_nlon-1)*(this->_nlat-1));
	e1v = NetCDF::readNetCDF(this->meshmask().c_str(),"e1v",1*(this->_nlon-1)*(this->_nlat-1));
	e1f = NetCDF::readNetCDF(this->meshmask().c_str(),"e1f",1*(this->_nlon-1)*(this->_nlat-1));

	double *e2t, *e2u, *e2v, *e2f;
	e2t = NetCDF::readNetCDF(this->meshmask().c_str(),"e2t",1*(this->_nlon-1)*(this->_nlat-1));
	e2u = NetCDF::readNetCDF(this->meshmask().c_str(),"e2u",1*(this->_nlon-1)*(this->_nlat-1));
	e2v = NetCDF::readNetCDF(this->meshmask().c_str(),"e2v",1*(this->_nlon-1)*(this->_nlat-1));
	e2f = NetCDF::readNetCDF(this->meshmask().c_str(),"e2f",1*(this->_nlon-1)*(this->_nlat-1));

	double *e3t, *e3u, *e3v, *e3w;
	e3t = NetCDF::readNetCDF(this->meshmask().c_str(),"e3t_0",this->_ncells);
	e3u = NetCDF::readNetCDF(this->meshmask().c_str(),"e3u_0",this->_ncells);
	e3v = NetCDF::readNetCDF(this->meshmask().c_str(),"e3v_0",this->_ncells);
	e3w = NetCDF::readNetCDF(this->meshmask().c_str(),"e3w_0",this->_ncells);

	// Allocate the fields
	this->_e1.set_dim(this->_ncells,4);
	this->_e2.set_dim(this->_ncells,4);
	this->_e3.set_dim(this->_ncells,4);

	// Fill the fields
	for (int kk = 0; kk < this->_nlev-1; kk++) {
		for (int jj = 0; jj < this->_nlat-1; jj++) {
			for (int ii = 0; ii < this->_nlon-1; ii++) {
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

	delete [] e1t; delete [] e1u; delete [] e1v; delete [] e1f;
	delete [] e2t; delete [] e2u; delete [] e2v; delete [] e2f;
	delete [] e3t; delete [] e3u; delete [] e3v; delete [] e3w;
}

std::string OGS::var_WritePath(int i, int j, const char *str, const char *token) {
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
	std::printf("    Mesh file:      %s\n", this->_meshfile.c_str());
	std::printf("    Meshmask file:  %s\n", this->_meshmask.c_str());
	std::printf("    Dimensions:     (%d,%d,%d)\n", this->_nlon,this->_nlat,this->_nlev);
	std::printf("    Longitude [m]:  (%.2f,%.2f)\n", this->lon2meters(0), this->lon2meters(-1));
	std::printf("    Latitude [m]:   (%.2f,%.2f)\n", this->lat2meters(0),this->lat2meters(-1));
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
	OGS *newOGS(const char *fname, const int nlon, const int nlat, const int nlev,
		double *lon2m, double *lat2m, double *nav_lev, double *bmask, double *cmask) {
		// Create a new instance of the class
		OGS *ogscls; ogscls = new OGS[1];
		// Populate the class
		ogscls->SetMfile(fname);
		ogscls->SetWdir(".");
		ogscls->Setlon2m(nlon,lon2m);
		ogscls->Setlat2m(nlat,lat2m);
		ogscls->Setnavlev(nlev,nav_lev);
		ogscls->Setncells();
		// Load the masks
		ogscls->SetMask(0,bmask);
		ogscls->SetMask(1,cmask);
		// Return
		return ogscls;
	}
	int OGSWriteMesh(OGS *ogscls) { return ogscls->writeMesh(); }
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