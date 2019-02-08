/*=========================================================================

  Module:    OGS File functions

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef OGSfile_h
#define OGSfile_h

#include "OGSdefs.h"

/* Function prototypes */

extern "C" int readOGSFile(const char *fname, char *mesh_file, char *mesh_mask,
	ave_var *ave_phys, ave_var *ave_freq, ogs_time *timeStepInfo);

extern "C" char *writeOGSPath(const char *str1, const char *str2, const char *token);

#endif