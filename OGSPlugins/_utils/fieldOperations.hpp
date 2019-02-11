/*=========================================================================

  Program:   Field Operations
  Module:    fieldOperations.hpp

  Useful operations with field arrays.

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef FIELDOPERATIONS_H
#define FIELDOPERATIONS_H

#include "field.h"

#define PNTIND(ii,jj,kk,nx,ny) ( (nx)*(ny)*(kk) + (nx)*(jj) + (ii) )

namespace field
{

	/* FACE2CELLPROJ

		Projection of a face centered array to a cell centered array. Useful
		for the velocity field or when computing gradients.

	*/
	template<class T>
	void Face2CellProj_ijk(int ii, int jj, int kk, int nx, int ny, 
		Field<T> &f, Field<T> &e1, Field <T> &e2, Field <T> &e3, T *out) {

		T w1, w2, w3;

		/* X COORDINATE */

		if (ii > 0) {
			w1 = e2[PNTIND(ii-1,jj,kk,nx,ny)][1]*e3[PNTIND(ii-1,jj,kk,nx,ny)][1]; // (e2u*e3u)_(i-1)
			w2 = e2[PNTIND(ii,jj,kk,nx,ny)][1]*e3[PNTIND(ii,jj,kk,nx,ny)][1];     // (e2u*e3u)_i
			w3 = e2[PNTIND(ii,jj,kk,nx,ny)][0]*e3[PNTIND(ii,jj,kk,nx,ny)][0];     // e2t*e3t

			out[0] = 0.5*( w1*f[PNTIND(ii-1,jj,kk,nx,ny)][0] + w2*f[PNTIND(ii,jj,kk,nx,ny)][0])/w3; // 0.5*(a1_(i-1)+a1_(i))/e2t/e3t
		}

		/* Y COORDINATE */

		if (jj > 0) {
			w1 = e1[PNTIND(ii-1,jj,kk,nx,ny)][1]*e3[PNTIND(ii-1,jj,kk,nx,ny)][1]; // (e1u*e3u)_(i-1)
			w2 = e1[PNTIND(ii,jj,kk,nx,ny)][1]*e3[PNTIND(ii,jj,kk,nx,ny)][1];     // (e1u*e3u)_i
			w3 = e1[PNTIND(ii,jj,kk,nx,ny)][0]*e3[PNTIND(ii,jj,kk,nx,ny)][0];     // e1t*e3t

			out[1] = 0.5*( w1*f[PNTIND(ii,jj-1,kk,nx,ny)][1] + w2*f[PNTIND(ii,jj,kk,nx,ny)][1])/w3; // 0.5*(a2_(j-1)+a2_(j))/e1t/e3t
		}

		/* Z COORDINATE */

		if (kk > 0)
			out[2] = 0.5*(f[PNTIND(ii,jj,kk-1,nx,ny)][2] + f[PNTIND(ii,jj,kk,nx,ny)][2]);
	}
	template<class T>
	void Face2CellProj(Field<T> &f, Field<T> e1, Field <T> e2, Field <T> e3, int nx, int ny, int nz) {

		T *out; out = new T[f.get_sz()];

		// Loop the components
		for (int kk = 0; kk < nz; kk++) {
			for (int jj = 0; jj < ny; jj++) { 
				for (int ii = 0; ii < nx; ii++) {
					Face2CellProj_ijk(ii,jj,kk,nx,ny,f,e1,e2,e3,
						out + f.get_m()*PNTIND(ii,jj,kk,nx,ny));
				} 
			}
		}

		f.set_val(out);
		delete [] out;
	}
}

#endif