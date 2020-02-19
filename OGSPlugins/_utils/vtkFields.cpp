/*=========================================================================

  Program:   VTK Fields
  Module:    vtkFields.cpp

  Wrapper to create many kinds of VTK fields 
  that might be needed accross plugins.

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include <cstdlib>
#include <cstring>

#ifdef __linux__
// Include OpenMP when working with GCC
#include <omp.h>
#define OMP_NUM_THREADS omp_get_num_threads()
#define OMP_THREAD_NUM  omp_get_thread_num()
#else
#define OMP_NUM_THREADS 1
#define OMP_THREAD_NUM  0
#endif

#include "vtkFields.h"

namespace VTK
{

	/* CREATEVTKSTRF

		Creates a vtk string array given the name 
		and data to be converted.

	*/
	vtkStringArray *createVTKstrf(const char *name, int n, const char *data) {
		// Create string array
		vtkStringArray *vtkArray = vtkStringArray::New();
		vtkArray->SetName(name);
		vtkArray->SetNumberOfTuples(n);

		// Set the value
		if (data != NULL) 
			for (int ii = 0; ii < n; ii++) vtkArray->SetValue(ii,data);

		return vtkArray;
	}

	/* CREATERECTILINEARGRID

		Creates a VTK rectilinear grid given the mesh dimensions and the 
		conversion to meters.

	*/
	void createRectilinearGrid(int nx, int ny, int nz, 
		double *x, double *y, double *z, double scalf, vtkRectilinearGrid *rgrid) {

		// Set dimension arrays
		vtkDoubleArray *vtkx; vtkx = createVTKscaf<vtkDoubleArray,double>("x coord", nx, x);
		vtkDoubleArray *vtky; vtky = createVTKscaf<vtkDoubleArray,double>("y coord", ny, y);
		vtkDoubleArray *vtkz; vtkz = createVTKscaf<vtkDoubleArray,double>("z coord", nz, z);

		// Fix scaling in z
		#pragma omp parallel
		{
		for (int ii=OMP_THREAD_NUM; ii<nz; ii+=OMP_NUM_THREADS) 
			vtkz->SetTuple1(ii,-scalf*z[ii]);
		}

		// Set rectilinear grid
		rgrid->SetDimensions(nx,ny,nz);
		rgrid->SetXCoordinates(vtkx); vtkx->Delete();
		rgrid->SetYCoordinates(vtky); vtky->Delete();
		rgrid->SetZCoordinates(vtkz); vtkz->Delete();
	}

}