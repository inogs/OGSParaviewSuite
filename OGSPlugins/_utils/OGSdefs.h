/*=========================================================================

  Module:    OGS Definitions

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef OGS_DEFS_h
#define OGS_DEFS_h

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
#define ERROR(errstr,errid) {fprintf(stderr,"Error in %s line %d:\n%s\n",__FILE__,__LINE__,errstr); return(-errid);}

/* ENUMERATIONS */

typedef struct _OGS_VAR { // OGS variable structure
	char name[VARSZ];	      // Variable name
  char vname[VARSZ];      // Variable name inside NetCDF
	char path[LBUFFSZ];	    // Path to file
}ogs_var;

typedef struct _AVE_VAR { // AVE_FREQ / AVE_PHYS structure
	int nvars = 0;          // Number of variables
	ogs_var *vars;		      // Variable structure
}ave_var;

typedef struct _OGS_TIME { // OGS timestep structure
  int ntsteps = 0;         // Number of time steps
  char **datetime;         // Date and hour
}ogs_time;

#endif