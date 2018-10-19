/*=========================================================================

  Module:    OGS File functions

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include <cstdlib>
#include <cstdio>
#include <string>
#include <cstring>

#include "OGSdefs.h"
#include "OGSfile.h"

using namespace std;

inline char *reads(char *line, int size, FILE *fin) {
	char *l = fgets(line,size,fin);
	if (l != NULL)
	  for (int ii = 0; line[ii] != '\0'; ii++)
	    if (line[ii] == 10) { line[ii] = '\0'; break; }
	return l;
}

inline char *trim(char *str) {
        char *end;
        while(isspace(*str)) str++; // Trim leading space
        if(*str == 0) return str;// All spaces?
        // Trim trailing space
        end = str + strlen(str) - 1;
        while(end > str && isspace(*end)) end--;
        // Write new null terminator
        *(end+1) = 0;
        return str;
}

/* OGS File reader
	TODO: description
*/
extern "C" int readOGSFile(const char *fname, char *mesh_file, char *mesh_mask,
	ave_var *ave_phys, ave_var *ave_freq, ogs_time *timeStepInfo) {
	
	char line[BUFFSZ], wrkdir[BUFFSZ];
	char *linep;
	int sec = 0, n_ave_phys_read = -1, n_ave_freq_read = -1, n_time_read = -1;
	int iwrkdir = 0, imesh = 0, iave_phys = 0, iave_freq = 0, itime = 0;

	// Open file for reading
	FILE *myfile;
	myfile = fopen(fname,"r"); if (myfile == NULL) ERROR("Cannot open file",1)

	// Read file line by line
	while(reads(line,sizeof(line),myfile)) {
		// Ingnore the lines with a comment
		if (line[0] == '#') continue;
		// Read keywords of sections
		if (string(line) == "WRKDIR")   {sec++; continue;}
		if (string(line) == "MESH")     {sec++; continue;}
		if (string(line) == "AVE_PHYS") {sec++; continue;}
		if (string(line) == "AVE_FREQ") {sec++; continue;}
		if (string(line) == "TIME")     {sec++; continue;}
		// Read the sections
		if (sec == 1 && !iwrkdir) { // Work directory
			strcpy(wrkdir,trim(line));
			iwrkdir = 1;
		}
		if (sec == 2 && !imesh) { // Mesh section, line contains path to the mesh file
			sprintf(mesh_file,"%s/%s",wrkdir,trim(line));
			sprintf(mesh_mask,"%s/meshmask.nc",wrkdir);
			imesh = 1;
		}
		if (sec == 3 && !iave_phys) { // Physical variables
			if (n_ave_phys_read >= 0) {
				if (n_ave_phys_read > ave_phys->nvars-1) {iave_phys = 1; continue;}
				// Split the string, read variable name
				strtok(line,":");
				strcpy(ave_phys->vars[n_ave_phys_read].name,trim(line));
				// Split the string, read netcfd name
				linep = strtok(NULL,":");
				strcpy(ave_phys->vars[n_ave_phys_read].vname,trim(linep));
				// Split the string, read variable path
				linep = strtok(NULL," ");
				sprintf(ave_phys->vars[n_ave_phys_read].path,"%s/%s",wrkdir,trim(linep));
				n_ave_phys_read++;
			} else {
				// Read number of variables
				ave_phys->nvars = atoi(line);
				// Allocate output array
				ave_phys->vars = (ogs_var*)malloc(ave_phys->nvars*sizeof(ogs_var));
				n_ave_phys_read++;
			}
		}
		if (sec == 4 && !iave_freq) { // Biogeochemical variables
			if (n_ave_freq_read >= 0) {
				if (n_ave_freq_read > ave_freq->nvars-1) {iave_freq = 1; continue;}
				// Split the string, read variable name
				strtok(line,":");
				strcpy(ave_freq->vars[n_ave_freq_read].name,trim(line));
				// Split the string, read netcfd name
				linep = strtok(NULL,":");
				strcpy(ave_freq->vars[n_ave_freq_read].vname,trim(linep));
				// Split the string, read variable path
				linep = strtok(NULL," ");
				sprintf(ave_freq->vars[n_ave_freq_read].path,"%s/%s",wrkdir,trim(linep));
				n_ave_freq_read++;
			} else {
				// Read number of variables
				ave_freq->nvars = atoi(line);
				// Allocate output array
				ave_freq->vars = (ogs_var*)malloc(ave_freq->nvars*sizeof(ogs_var));
				n_ave_freq_read++;
			}
		}
		if (sec == 5 && !itime) { // Timestep information
			if (n_time_read >= 0) {
				if (n_time_read > timeStepInfo->ntsteps-1) {itime = 1; continue;}
				// Read datetime
				strcpy(timeStepInfo->datetime[n_time_read],trim(line));
				n_time_read++;
			} else {
				// Read number of timesteps
				timeStepInfo->ntsteps = atoi(line);
				// Allocate
				timeStepInfo->datetime = (char**)malloc(timeStepInfo->ntsteps*sizeof(char*));
				for (int ii = 0; ii < timeStepInfo->ntsteps; ii++)
					timeStepInfo->datetime[ii] = (char*)malloc(VARSZ*CHRSZ);
				n_time_read++;
			}

		}
	}
	// Close file
	fclose(myfile);

	return 1;
}

extern "C" char *writeOGSPath(const char *str1, const char *str2, const char *token) {
	char aux[VARSZ];
	char *fname;

	fname = (char*)malloc(LBUFFSZ*CHRSZ);

	strcpy(aux,str2);
	strcpy(fname,str1);

	// Split by token and concatenate str2
	strtok(fname,token);
	strcat(aux,trim(strtok(NULL,token)));

	// Build path
	strcat(fname,aux);

	return(fname);
}
