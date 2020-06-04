/*=========================================================================

  Module:    Macros

  Define macros for all the ParaView filters. 

  Copyright (c) 2020 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkFloatArray.h"
#include "vtkDoubleArray.h"
#include "vtkTypeUInt8Array.h"
#include <cstdint>

/* ARRAY PRECISION

	The macros here define the array precision for all the
	ParaView filters.
		> FLDARRAY: Field arrays precision
		> VTKARRAY: VTK array precision
		> FLDMASK:  Mask array precision (uint8_t = 8 byte)
		> VTKMASK:  VTK mask array precision
*/

#ifdef VTK_USE_DOUBLE
#define FLDARRAY double
#define VTKARRAY vtkDoubleArray
#else
#define FLDARRAY float
#define VTKARRAY vtkFloatArray
#endif

#define FLDMASK  uint8_t
#define VTKMASK  vtkTypeUInt8Array

/* OPENMP SUPPORT

  Sets up macros for OpenMP parallelization.
*/

#ifdef USE_OMP
#include <omp.h>
#define OMP_THREAD_NUM  omp_get_thread_num()
#define OMP_NUM_THREADS omp_get_num_threads()
#define OMP_MAX_THREADS omp_get_max_threads()
#else
#define OMP_THREAD_NUM  0
#define OMP_NUM_THREADS 1
#define OMP_MAX_THREADS 1
#endif