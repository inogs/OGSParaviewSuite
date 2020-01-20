/*=========================================================================

  Program:   VTK Operations
  Module:    vtkOperations.cpp

  Wrapper to perform operations with VTK 
  that might be needed accross plugins.

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

	 This software is distributed WITHOUT ANY WARRANTY; without even
	 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
	 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkCell.h"
#include "vtkGenericCell.h"
#include "vtkPoints.h"
#include "vtkCellData.h"
#include "vtkPointData.h"
#include "vtkFieldData.h"
#include "vtkSmartPointer.h"

#ifdef __linux__
// Include OpenMP when working with GCC
#include <omp.h>
#define OMP_NUM_THREADS omp_get_num_threads()
#define OMP_THREAD_NUM  omp_get_thread_num()
#else
#define OMP_NUM_THREADS 1
#define OMP_THREAD_NUM  0
#endif

#include "vtkOperations.h"

#define CLLIND(ii,jj,kk,nx,ny) ( (nx-1)*(ny-1)*(kk) + (nx-1)*(jj) + (ii) )
#define PNTIND(ii,jj,kk,nx,ny) ( (nx)*(ny)*(kk) + (nx)*(jj) + (ii) )

namespace VTK
{
	/* GETVTKCELLPOINTS
	
		Recovers the cell point coordinates given a VTKDataSet.

	*/
	v3::V3v getVTKCellPoints(vtkDataSet *mesh, double fact) {
		// Define output array
		v3::V3v xyz( mesh->GetNumberOfPoints() );

		// Loop the mesh and get the point coordinates
		#pragma omp parallel shared(mesh,xyz)
		{
		for (int ii=OMP_THREAD_NUM; ii < xyz.len(); ii+=OMP_NUM_THREADS) {
			// Obtain point from VTK structure
			double pnt[3]; mesh->GetPoint(ii,pnt);
			// Set the V3 structure
			xyz[ii][0] = pnt[0];
			xyz[ii][1] = pnt[1];
			xyz[ii][2] = pnt[2]/fact;
		}
		}

		// Return V3v
		return xyz;
	}

	/* GETVTKCELLCENTERS

		Recovers the cell center coordinates given a VTKDataSet.

	*/
	v3::V3v getVTKCellCenters(vtkDataSet *mesh, double fact) {
		// Define output array
		v3::V3v xyz( mesh->GetNumberOfCells() );

		// Loop the mesh and get the point coordinates
		#pragma omp parallel shared(mesh,xyz) firstprivate(fact)
		{
		vtkSmartPointer<vtkGenericCell> cell = vtkSmartPointer<vtkGenericCell>::New();
		for (int ii=OMP_THREAD_NUM; ii < xyz.len(); ii+=OMP_NUM_THREADS) {
			// Preallocate to zero
			xyz[ii] = v3::V3(0.,0.,0.);
			// Get number of points per cell
			mesh->GetCell(ii,cell);
			int ncellp    = cell->GetNumberOfPoints();
			// Loop the points in the cell
			for (int pId = 0; pId < ncellp; ++pId) {
				// Obtain point from VTK structure
				double pnt[3]; cell->GetPoints()->GetPoint(pId,pnt);
				// Set the V3 structure
				xyz[ii][0] += pnt[0]/ncellp;
				xyz[ii][1] += pnt[1]/ncellp;
				xyz[ii][2] += pnt[2]/ncellp;						
			}
			xyz[ii][2] /= fact;
		}
		}

		// Return V3v
		return xyz;		
	}

}