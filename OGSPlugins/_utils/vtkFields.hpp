/*=========================================================================

  Program:   VTK Fields
  Module:    vtkFields.hpp

  Wrapper to create many kinds of VTK fields 
  that might be needed accross plugins.

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef VTKFIELDS_H
#define VTKFIELDS_H

#include "vtkRectilinearGrid.h"
#include "vtkFloatArray.h"
#include "vtkDoubleArray.h"
#include "vtkStringArray.h"

#include "field.h"

#define CLLIND(ii,jj,kk,nx,ny)          ( (nx-1)*(ny-1)*(kk) + (nx-1)*(jj) + (ii) )
#define PNTIND(ii,jj,kk,nx,ny)          ( (nx)*(ny)*(kk) + (nx)*(jj) + (ii) )

namespace VTK
{

	/* CREATEVTKFROMFIELD
		
		Creates a vtk array from a field given the array name.
	*/
	template <class ARRAY, class P>
	ARRAY *createVTKfromField(const char *name, field::Field<P> &f) {
		// Create VTK array
		ARRAY *vtkArray = ARRAY::New();
		vtkArray->SetName(name);
		vtkArray->SetNumberOfComponents(f.get_m()); // Scalar field
		vtkArray->SetNumberOfTuples(f.get_n());

		if (f.data() != NULL)
			// Fill the vtkArray with the values of the array
			std::memcpy(vtkArray->GetPointer(0),f.data(),f.get_sz()*sizeof(P));
		else
			// Preallocate vtkArray to zero
			vtkArray->Fill(0.);

		return vtkArray;
	}

	/* CREATEFIELDFROMVTK
		
		Creates a field froma vtk array.
	*/
	template <class ARRAY, class P>
	field::Field<P> createFieldfromVTK(ARRAY *vtkArray) {
		// Recover dimensions from VTK
		int n = vtkArray->GetNumberOfTuples();
		int m = vtkArray->GetNumberOfComponents();
		// Create field
		field::Field<P> f(n,m,vtkArray->GetPointer(0));

		return f;
	}

	/* CREATEVTKSCAF
		
		Creates an array for a scalar field given the array name, its dimensions
		and an additional array to fill.
	*/
	template <class ARRAY, class P>
	ARRAY *createVTKscaf(const char *name, int n, P *array) {
		// Create VTK array
		ARRAY *vtkArray = ARRAY::New();
		vtkArray->SetName(name);
		vtkArray->SetNumberOfComponents(1); // Scalar field
		vtkArray->SetNumberOfTuples(n);
		
		if (array != NULL)
			// Fill the vtkArray with the values of the array
			std::memcpy(vtkArray->GetPointer(0),array,n*sizeof(P));
		else
			// Preallocate vtkArray to zero
			vtkArray->Fill(0.);

		return vtkArray;
	}
	template <class ARRAY, class P>
	ARRAY *createVTKscaf(const char *name, int nx, int ny, int nz, P *array) {
		// Create float array
		ARRAY *vtkArray = ARRAY::New();
		vtkArray->SetName(name);
		vtkArray->SetNumberOfComponents(1); // Scalar field
		vtkArray->SetNumberOfTuples(nx*ny*nz);
		
		if (array != NULL) {
			// Fill the vtkArray with the values of the array
			for (int kk = 0; kk < nz; kk++) {
				for (int jj = 0; jj < ny; jj++) {
					for (int ii = 0; ii < nx; ii++) {
						vtkArray->SetTuple1(PNTIND(ii,jj,kk,nx,ny),
							array[PNTIND(ii,jj,kk,nx,ny)]);
					}
				}
			}
		} else {
			// Preallocate vtkArray to zero
			vtkArray->Fill(0.); 
		}

		return vtkArray;
	}

	/* CREATEVTKSTRF

		Creates a vtk string array given the name 
		and data to be converted.

	*/
	vtkStringArray *createVTKstrf(const char *name, int n, const char *data);

	/* CREATERECTILINEARGRID

		Creates a VTK rectilinear grid given the mesh dimensions and the 
		conversion to meters.

	*/
	void createRectilinearGrid(int nx, int ny, int nz, 
		double *x, double *y, double *z, double scalf, vtkRectilinearGrid *rgrid);

}

#endif