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
#include "vtkDataSet.h"
#include "vtkSmartPointer.h"

#ifdef USE_OMP
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

	/* GETCONNECTEDPOINTS

		Recovers all the points connected to a point in the mesh.
	*/
	std::vector<int> getConnectedPoints(vtkDataSet *mesh, const int pointId) {
		
		std::vector<int> connectedPoints;
		vtkIdList *cellIdList = vtkIdList::New();

		// Get all the cells that belong to the current point
		mesh->GetPointCells(pointId,cellIdList);

		// Loop all the cells that we just found
		for (int cellId = 0; cellId < cellIdList->GetNumberOfIds(); ++cellId) {
			vtkIdList *pointIdList = vtkIdList::New();
			// Recover points in the cell
			mesh->GetCellPoints(cellIdList->GetId(cellId),pointIdList);
			// Append to connected
			if (pointIdList->GetId(0) != pointId)
				connectedPoints.push_back( pointIdList->GetId(0) );
			else
				connectedPoints.push_back( pointIdList->GetId(1) );
			// Delete
			pointIdList->Delete();
		}

		cellIdList->Delete(); 
		return connectedPoints;
	}

	/* GETCONNECTEDCELLS

		Recovers all the cells connected to a cell in the mesh.
	*/
	std::vector<int> getConnectedCells(vtkDataSet *mesh, const int cellId) {

		std::vector<int> connectedCells;
		vtkIdList *pointIdList = vtkIdList::New();
		mesh->GetCellPoints(cellId,pointIdList);

		for (int pointId = 0; pointId<pointIdList->GetNumberOfIds(); ++pointId) {
			// Create a list with the points of the edge of the cell
			vtkIdList *edgeIdList = vtkIdList::New();
			// Add one edge
			edgeIdList->InsertNextId(pointIdList->GetId(pointId));
			// Add the other edge
			if (pointId+1 == pointIdList->GetNumberOfIds())
				edgeIdList->InsertNextId(pointIdList->GetId(0));
			else
				edgeIdList->InsertNextId(pointIdList->GetId(pointId+1));
			// Get the neighbors of the cell
			vtkIdList *cellIdList = vtkIdList::New();
			mesh->GetCellNeighbors(cellId,edgeIdList,cellIdList);
			// Load into the vector
			for (int cId = 0; cId<cellIdList->GetNumberOfIds(); ++cId)
				connectedCells.push_back( cellIdList->GetId(cId) );
			// Delete
			edgeIdList->Delete(); cellIdList->Delete();
		}

		pointIdList->Delete();
		return connectedCells;
	}

}