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
#include "vtkPoints.h"
#include "vtkCellData.h"
#include "vtkPointData.h"
#include "vtkFieldData.h"

#include "vtkOperations.hpp"

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
		v3::V3v::iterator iter;
		for(iter = xyz.begin(); iter != xyz.end(); ++iter) {
			// Obtain point from VTK structure
			double pnt[3]; mesh->GetPoint(iter.ind(),pnt);
			// Set the V3 structure
			iter[0] = pnt[0];
			iter[1] = pnt[1];
			iter[2] = pnt[2]/fact;
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
		v3::V3v::iterator iter;
		for(iter = xyz.begin(); iter != xyz.end(); ++iter) {
			// Get number of points per cell
			vtkCell *cell = mesh->GetCell(iter.ind());
			int ncellp    = cell->GetNumberOfPoints();
			// Preallocate to zero
			*iter = v3::V3(0.,0.,0.);
			// Loop the points in the cell
			for (int pId = 0; pId < ncellp; ++pId) {
				// Obtain point from VTK structure
				double pnt[3]; cell->GetPoints()->GetPoint(pId,pnt);
				// Set the V3 structure
				iter[0] += pnt[0]/ncellp;
				iter[1] += pnt[1]/ncellp;
				iter[2] += pnt[2]/ncellp;						
			}
			iter[2] /= fact;

		}

		// Return V3v
		return xyz;		
	}

}