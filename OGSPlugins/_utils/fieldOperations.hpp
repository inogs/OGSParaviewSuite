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

#include "V3.h"
#include "field.h"

#include <vector>

#define PNTIND(ii,jj,kk,nx,ny) ( (nx)*(ny)*(kk) + (nx)*(jj) + (ii) )

namespace field
{

	/* UVW2T

		Projection of a face centered array to a cell centered array. Useful
		for the velocity field or when computing gradients.

	*/
	template<class T>
	void UVW2T_ijk(int ii, int jj, int kk, int nx, int ny, T f1[], T f2[],
		Field<T> &e1, Field <T> &e2, Field <T> &e3, T *out) {

		T w1, w2, w3;

		/* X COORDINATE */

		if (ii > 0) {
			w1 = e2[PNTIND(ii-1,jj,kk,nx,ny)][1]*e3[PNTIND(ii-1,jj,kk,nx,ny)][1]; // (e2u*e3u)_(i-1)
			w2 = e2[PNTIND(ii,jj,kk,nx,ny)][1]*e3[PNTIND(ii,jj,kk,nx,ny)][1];     // (e2u*e3u)_i
			w3 = e2[PNTIND(ii,jj,kk,nx,ny)][0]*e3[PNTIND(ii,jj,kk,nx,ny)][0];     // e2t*e3t

			out[0] = 0.5*( w1*f1[0] + w2*f2[0])/w3; // 0.5*(a1_(i-1)+a1_(i))/e2t/e3t
		} else {
			out[0] = f2[0];
		}

		/* Y COORDINATE */

		if (jj > 0) {
			w1 = e1[PNTIND(ii-1,jj,kk,nx,ny)][1]*e3[PNTIND(ii-1,jj,kk,nx,ny)][1]; // (e1u*e3u)_(i-1)
			w2 = e1[PNTIND(ii,jj,kk,nx,ny)][1]*e3[PNTIND(ii,jj,kk,nx,ny)][1];     // (e1u*e3u)_i
			w3 = e1[PNTIND(ii,jj,kk,nx,ny)][0]*e3[PNTIND(ii,jj,kk,nx,ny)][0];     // e1t*e3t

			out[1] = 0.5*( w1*f1[1] + w2*f2[1])/w3; // 0.5*(a2_(j-1)+a2_(j))/e1t/e3t
		} else {
			out[1] = f2[1];
		}

		/* Z COORDINATE */

		if (kk > 0)
			out[2] = 0.5*(f1[2] + f2[2]);
		else
			out[2] = f2[2];
	}
	template<class T>
	Field<T> UVW2T(Field<T> &f_UVW, Field<T> &e1, Field <T> &e2, Field <T> &e3, int nx, int ny, int nz) {
		// Create new field
		Field<T> f_T(f_UVW);
		// Loop the components
		for (int kk = 0; kk < nz; ++kk) {
			for (int jj = 0; jj < ny; ++jj) { 
				for (int ii = 0; ii < nx; ++ii) {
					int ind  = PNTIND(ii,jj,kk,nx,ny);
					int indx = PNTIND(ii-1,jj,kk,nx,ny);
					int indy = PNTIND(ii,jj-1,kk,nx,ny);
					int indz = PNTIND(ii,jj,kk-1,nx,ny);

					T f1_UVW[3] = {f_UVW[indx][0],f_UVW[indy][1],f_UVW[indz][2]};
					T f2_UVW[3] = {f_UVW[ind][0], f_UVW[ind][1], f_UVW[ind][2]};

					UVW2T_ijk(ii,jj,kk,nx,ny,f1_UVW,f2_UVW,e1,e2,e3,f_T[ind]);
				} 
			}
		}
		return f_T;
	}

	/* GRADXYZ2

		Second order, face centered, approximation of the gradient
		using a standard, finite difference approach.

	*/
	template<class T>
	void gradXYZ2_ijk(int ii, int jj, int kk, int nx, int ny, int nz,
		v3::V3v &xyz, Field<T> &f, T *deri) {

		int ind, ind1;

		/* DERIVATIVES WITH RESPECT TO X */
		if (ii == 0) {
			ind  = PNTIND(ii+1,jj,kk,nx,ny); ind1 = PNTIND(ii,jj,kk,nx,ny);
			ind  = PNTIND(ii+1,jj,kk,nx,ny); ind1 = PNTIND(ii,jj,kk,nx,ny);
		} else if (ii == nx-1) {
			ind  = PNTIND(ii,jj,kk,nx,ny);   ind1 = PNTIND(ii-1,jj,kk,nx,ny);
			ind  = PNTIND(ii,jj,kk,nx,ny);   ind1 = PNTIND(ii-1,jj,kk,nx,ny);
		} else {
			ind  = PNTIND(ii+1,jj,kk,nx,ny); ind1 = PNTIND(ii-1,jj,kk,nx,ny);
			ind  = PNTIND(ii+1,jj,kk,nx,ny); ind1 = PNTIND(ii-1,jj,kk,nx,ny);
		}
		// Compute the derivatives with respect to x
		for (int ii = 0; ii < f.get_m(); ++ii) {
			deri[0 + 3*ii] = (f[ind][ii] - f[ind1][ii])/(xyz[ind][0] - xyz[ind1][0]); // dqdx
		}

		/* DERIVATIVES WITH RESPECT TO Y */
		if (jj == 0) {
			ind  = PNTIND(ii,jj+1,kk,nx,ny); ind1 = PNTIND(ii,jj,kk,nx,ny);
			ind  = PNTIND(ii,jj+1,kk,nx,ny); ind1 = PNTIND(ii,jj,kk,nx,ny);
		} else if (jj == ny-1) {
			ind  = PNTIND(ii,jj,kk,nx,ny);   ind1 = PNTIND(ii,jj-1,kk,nx,ny);
			ind  = PNTIND(ii,jj,kk,nx,ny);   ind1 = PNTIND(ii,jj-1,kk,nx,ny);
		} else {
			ind  = PNTIND(ii,jj+1,kk,nx,ny); ind1 = PNTIND(ii,jj-1,kk,nx,ny);
			ind  = PNTIND(ii,jj+1,kk,nx,ny); ind1 = PNTIND(ii,jj-1,kk,nx,ny);
		}
		// Compute the derivatives with respect to x
		for (int ii = 0; ii < f.get_m(); ++ii)
			deri[1 + 3*ii] = (f[ind][ii] - f[ind1][ii])/(xyz[ind][1] - xyz[ind1][1]); // dqdy

		/* DERIVATIVES WITH RESPECT TO Z */
		if (kk == 0) {
			ind  = PNTIND(ii,jj,kk+1,nx,ny); ind1 = PNTIND(ii,jj,kk,nx,ny);
			ind  = PNTIND(ii,jj,kk+1,nx,ny); ind1 = PNTIND(ii,jj,kk,nx,ny);
		} else if (kk == nz-1) {
			ind  = PNTIND(ii,jj,kk,nx,ny);   ind1 = PNTIND(ii,jj,kk-1,nx,ny);
			ind  = PNTIND(ii,jj,kk,nx,ny);   ind1 = PNTIND(ii,jj,kk-1,nx,ny);
		} else {
			ind  = PNTIND(ii,jj,kk+1,nx,ny); ind1 = PNTIND(ii,jj,kk-1,nx,ny);
			ind  = PNTIND(ii,jj,kk+1,nx,ny); ind1 = PNTIND(ii,jj,kk-1,nx,ny);
		}
		// Compute the derivatives with respect to x
		for (int ii = 0; ii < f.get_m(); ++ii)
			deri[2 + 3*ii] = (f[ind][ii] - f[ind1][ii])/(xyz[ind][2] - xyz[ind1][2]); // dqdz
	}
	template<class T>
	Field<T> gradXYZ2(int nx, int ny, int nz, v3::V3v &xyz, 
		Field<T> &f, Field<T> &div, Field<T> &curl, Field<T> &Q) {
		// Create output array
		Field<T> grad(f.get_n(),3*f.get_m());
		// Loop the components
		for (int kk = 0; kk < nz; ++kk) {
			for (int jj = 0; jj < ny; ++jj) { 
				for (int ii = 0; ii < nx; ++ii) {
					// Point id
					int ind = PNTIND(ii,jj,kk,nx,ny);
					// Compute the gradient
					gradXYZ2_ijk(ii,jj,kk,nx,ny,nz,xyz,f,grad[ind]);
					// Computation of the divergence
					if (!div.isempty()) {
						div[ind][0] = grad[ind][0] + grad[ind][4] + grad[ind][8];
					}
					// Computation of the curl
					if (!curl.isempty()) {
						curl[ind][0] = grad[ind][7] - grad[ind][5];
						curl[ind][1] = grad[ind][2] - grad[ind][6];
						curl[ind][2] = grad[ind][3] - grad[ind][1];
					}
					// Computation of the Q-criterion
					if (!Q.isempty()) {
						Q[ind][0] = -0.5*(grad[ind][0]*grad[ind][0] + grad[ind][4]*grad[ind][4] + grad[ind][8]*grad[ind][8])
						            -grad[ind][1]*grad[ind][3]-grad[ind][2]*grad[ind][6]-grad[ind][5]*grad[ind][7];
					}
				}
			}
		}
		return grad;
	}

	/* GRADXY4

		Fourth order, face centered, approximation of the gradient using
		a standard finite difference approach.

	*/
	template<class T>
	void gradXYZ4_ijk(int ii, int jj, int kk, int nx, int ny, int nz,
		v3::V3v &xyz, Field<T> &f, T *deri) {

		int ind, ind1, ind2, ind3;

		/* DERIVATIVES WITH RESPECT TO X */
		if (ii <= 1) {
			ind  = PNTIND(ii+2,jj,kk,nx,ny);
			ind1 = PNTIND(ii+1,jj,kk,nx,ny);
			ind2 = PNTIND(ii,jj,kk,nx,ny);

			for (int ii = 0; ii < f.get_m(); ++ii)
				deri[0 + 3*ii] = (-f[ind][ii] + 4.*f[ind1][ii] - 3.*f[ind2][ii])/2./(xyz[ind1][0] - xyz[ind2][0]);
		} else if (ii >= nx-2) {
			ind  = PNTIND(ii,jj,kk,nx,ny);
			ind1 = PNTIND(ii-1,jj,kk,nx,ny);
			ind2 = PNTIND(ii-2,jj,kk,nx,ny);

			for (int ii = 0; ii < f.get_m(); ++ii)
				deri[0 + 3*ii] = (3.*f[ind][ii] - 4.*f[ind1][ii] + f[ind2][ii])/2./(xyz[ind][0] - xyz[ind1][0]);
		} else {
			ind  = PNTIND(ii+2,jj,kk,nx,ny);
			ind1 = PNTIND(ii+1,jj,kk,nx,ny);
			ind2 = PNTIND(ii-1,jj,kk,nx,ny);
			ind3 = PNTIND(ii-2,jj,kk,nx,ny);

			for (int ii = 0; ii < f.get_m(); ++ii)
				deri[0 + 3*ii] = (-f[ind][ii] + 8.*f[ind1][ii] - 8.*f[ind2][ii] + f[ind3][ii])/6./(xyz[ind1][0] - xyz[ind2][0]);
		}
					
		/* DERIVATIVES WITH RESPECT TO Y */
		if (jj <= 1) {
			ind  = PNTIND(ii,jj+2,kk,nx,ny);
			ind1 = PNTIND(ii,jj+1,kk,nx,ny);
			ind2 = PNTIND(ii,jj,kk,nx,ny);

			for (int ii = 0; ii < f.get_m(); ++ii)
				deri[1 + 3*ii] = (-f[ind][ii] + 4.*f[ind1][ii] - 3.*f[ind2][ii])/2./(xyz[ind1][1] - xyz[ind2][1]);
		} else if (jj >= ny-2) {
			ind  = PNTIND(ii,jj,kk,nx,ny);
			ind1 = PNTIND(ii,jj-1,kk,nx,ny);
			ind2 = PNTIND(ii,jj-2,kk,nx,ny);

			for (int ii = 0; ii < f.get_m(); ++ii)
				deri[1 + 3*ii] = (3.*f[ind][ii] - 4.*f[ind1][ii] + f[ind2][ii])/2./(xyz[ind][1] - xyz[ind1][1]);
		} else {
			ind  = PNTIND(ii,jj+2,kk,nx,ny);
			ind1 = PNTIND(ii,jj+1,kk,nx,ny);
			ind2 = PNTIND(ii,jj-1,kk,nx,ny);
			ind3 = PNTIND(ii,jj-2,kk,nx,ny);

			for (int ii = 0; ii < f.get_m(); ++ii)
				deri[1 + 3*ii] = (-f[ind][ii] + 8.*f[ind1][ii] - 8.*f[ind2][ii] + f[ind3][ii])/6./(xyz[ind1][1] - xyz[ind2][1]);
		}
	
		/* DERIVATIVES WITH RESPECT TO Z */
		if (kk <= 1) {
			ind  = PNTIND(ii,jj,kk+2,nx,ny);
			ind1 = PNTIND(ii,jj,kk+1,nx,ny);
			ind2 = PNTIND(ii,jj,kk,nx,ny);

			for (int ii = 0; ii < f.get_m(); ++ii)
				deri[2 + 3*ii] = (-f[ind][ii] + 4.*f[ind1][ii] - 3.*f[ind2][ii])/2./(xyz[ind1][2] - xyz[ind2][2]);
		} else if (kk >= nz-2) {
			ind  = PNTIND(ii,jj,kk,nx,ny);
			ind1 = PNTIND(ii,jj,kk-1,nx,ny);
			ind2 = PNTIND(ii,jj,kk-2,nx,ny);

			for (int ii = 0; ii < f.get_m(); ++ii)
				deri[2 + 3*ii] = (3.*f[ind][ii] - 4.*f[ind1][ii] + f[ind2][ii])/2./(xyz[ind][2] - xyz[ind1][2]);
		} else {
			ind  = PNTIND(ii,jj,kk+2,nx,ny);
			ind1 = PNTIND(ii,jj,kk+1,nx,ny);
			ind2 = PNTIND(ii,jj,kk-1,nx,ny);
			ind3 = PNTIND(ii,jj,kk-2,nx,ny);

			for (int ii = 0; ii < f.get_m(); ++ii)
				deri[2 + 3*ii] = (-f[ind][ii] + 8.*f[ind1][ii] - 8.*f[ind2][ii] + f[ind3][ii])/6./(xyz[ind1][2] - xyz[ind2][2]);			
		}
	}
	template<class T>
	Field<T> gradXYZ4(int nx, int ny, int nz, v3::V3v &xyz, 
		Field<T> &f, Field<T> &div, Field<T> &curl, Field<T> &Q) {
		// Create output array
		Field<T> grad(f.get_n(),3*f.get_m());
		// Loop the components
		for (int kk = 0; kk < nz; ++kk) {
			for (int jj = 0; jj < ny; ++jj) { 
				for (int ii = 0; ii < nx; ++ii) {
					// Point id
					int ind = PNTIND(ii,jj,kk,nx,ny);
					// Compute the gradient
					gradXYZ4_ijk(ii,jj,kk,nx,ny,nz,xyz,f,grad[ind]);
					// Computation of the divergence
					if (!div.isempty()) {
						div[ind][0] = grad[ind][0] + grad[ind][4] + grad[ind][8];
					}
					// Computation of the curl
					if (!curl.isempty()) {
						curl[ind][0] = grad[ind][7] - grad[ind][5];
						curl[ind][1] = grad[ind][2] - grad[ind][6];
						curl[ind][2] = grad[ind][3] - grad[ind][1];
					}
					// Computation of the Q-criterion
					if (!Q.isempty()) {
						Q[ind][0] = -0.5*(grad[ind][0]*grad[ind][0] + grad[ind][4]*grad[ind][4] + grad[ind][8]*grad[ind][8])
						            -grad[ind][1]*grad[ind][3]-grad[ind][2]*grad[ind][6]-grad[ind][5]*grad[ind][7];
					}
				}
			}
		}
		return grad;
	}

	/* GRADOGS1

		First order, forward Euler, approximation of the gradient for
		a field using OGSTM-BFM approach. This gradient lives in the staggered 
		mesh, and it is projected back to the cell centered mesh. 

		The input is assumed to be on the cell centered mesh. 

	*/
	template<class T>
	void gradOGS1_ijk(int ii, int jj, int kk, int nx, int ny, int nz,
		Field<T> &f, T e1u, T e2v, T e3w, T *deri) {

		int ind, ind1;

		/* DERIVATIVES WITH RESPECT TO X */
		if (ii == nx-1) {
			ind  = PNTIND(ii,jj,kk,nx,ny); 
			ind1 = PNTIND(ii-1,jj,kk,nx,ny);
		} else {
			ind  = PNTIND(ii+1,jj,kk,nx,ny); 
			ind1 = PNTIND(ii,jj,kk,nx,ny);
		}
		for (int ii = 0; ii < f.get_m(); ++ii)
			deri[0 + 3*ii] = (f[ind][ii] - f[ind1][ii])/e1u; // dqdx

		/* DERIVATIVES WITH RESPECT TO Y */
		if (jj == ny-1) {
			ind  = PNTIND(ii,jj,kk,nx,ny); 
			ind1 = PNTIND(ii,jj-1,kk,nx,ny);
		} else {
			ind  = PNTIND(ii,jj+1,kk,nx,ny); 
			ind1 = PNTIND(ii,jj,kk,nx,ny);
		}
		for (int ii = 0; ii < f.get_m(); ++ii)
			deri[1 + 3*ii] = (f[ind][ii] - f[ind1][ii])/e2v; // dqdy

		/* DERIVATIVES WITH RESPECT TO Z */
		if (kk == nz-1) {
			ind  = PNTIND(ii,jj,kk,nx,ny); 
			ind1 = PNTIND(ii,jj,kk-1,nx,ny);
		} else {
			ind  = PNTIND(ii,jj,kk+1,nx,ny); 
			ind1 = PNTIND(ii,jj,kk,nx,ny);
		}
		for (int ii = 0; ii < f.get_m(); ++ii)
			deri[2 + 3*ii] = (f[ind][ii] - f[ind1][ii])/e3w; // dqdy
	}
	template<class T>
	Field<T> gradOGS1(int nx, int ny, int nz, Field<T> &f, 
		Field<T> &e1, Field<T> &e2, Field<T> &e3, Field<T> &div, Field<T> &curl, Field<T> &Q) {
		// Create output array
		Field<T> grad_UVW(f.get_n(),3*f.get_m(),0.), grad_T(f.get_n(),3*f.get_m());
		// Loop the components
		for (int kk = 0; kk < nz; ++kk) {
			for (int jj = 0; jj < ny; ++jj) { 
				for (int ii = 0; ii < nx; ++ii) {
					// Point id
					int ind = PNTIND(ii,jj,kk,nx,ny);
					int indx = PNTIND(ii-1,jj,kk,nx,ny);
					int indy = PNTIND(ii,jj-1,kk,nx,ny);
					int indz = PNTIND(ii,jj,kk-1,nx,ny);
					// Compute the gradient on UVW
					gradOGS1_ijk(ii,jj,kk,nx,ny,nz,f,e1[ind][1],e2[ind][2],e3[ind][3],grad_UVW[ind]);
					// Project from the UVW grid to T grid
					for (int gg = 0; gg < 3; ++gg) {
						T g1_UVW[3] = {grad_UVW[indx][0 + 3*gg],grad_UVW[indy][1 + 3*gg],grad_UVW[indz][2 + 3*gg]};
						T g2_UVW[3] = {grad_UVW[ind][0 + 3*gg], grad_UVW[ind][1 + 3*gg], grad_UVW[ind][2 + 3*gg]};
						UVW2T_ijk(ii,jj,kk,nx,ny,g1_UVW,g2_UVW,e1,e2,e3,grad_T[ind] + 3*gg);
					}
					// Computation of the divergence
					if (!div.isempty()) {
						div[ind][0] = grad_T[ind][0] + grad_T[ind][4] + grad_T[ind][8];
					}
					// Computation of the curl
					if (!curl.isempty()) {
						curl[ind][0] = grad_T[ind][7] - grad_T[ind][5];
						curl[ind][1] = grad_T[ind][2] - grad_T[ind][6];
						curl[ind][2] = grad_T[ind][3] - grad_T[ind][1];
					}
					// Computation of the Q-criterion
					if (!Q.isempty()) {
						Q[ind][0] = -0.5*(grad_T[ind][0]*grad_T[ind][0] + grad_T[ind][4]*grad_T[ind][4] + grad_T[ind][8]*grad_T[ind][8])
						            -grad_T[ind][1]*grad_T[ind][3]-grad_T[ind][2]*grad_T[ind][6]-grad_T[ind][5]*grad_T[ind][7];
					}
				}
			}
		}
		return grad_T;
	}

	/* GRADOGS2

		Second order, face centered, approximation of the gradient for
		a field using OGSTM-BFM approach. This gradient lives in the staggered 
		mesh, and it is projected back to the cell centered mesh. 

		The input is assumed to be on the cell centered mesh.
		This gradient method is experimental.

	*/
	template<class T>
	void gradOGS2_ijk(int ii, int jj, int kk, int nx, int ny, int nz,
		Field<T> &f, T e1u, T e2v, T e3w, T *deri) {

		int ind = 0, ind1 = 0;

		/* DERIVATIVES WITH RESPECT TO X */
		if (ii == 0) {
			ind  = PNTIND(ii+1,jj,kk,nx,ny); ind1 = PNTIND(ii,jj,kk,nx,ny);
		} else if (ii == nx-1) {
			ind  = PNTIND(ii,jj,kk,nx,ny);   ind1 = PNTIND(ii-1,jj,kk,nx,ny);
		} else {
			ind  = PNTIND(ii+1,jj,kk,nx,ny); ind1 = PNTIND(ii-1,jj,kk,nx,ny);
			e1u *= 2.;
		}
		for (int ii = 0; ii < f.get_m(); ++ii)
			deri[0 + 3*ii] = (f[ind][ii] - f[ind1][ii])/e1u; // dqdx

		/* DERIVATIVES WITH RESPECT TO Y */
		if (jj == 0) {
			ind  = PNTIND(ii,jj+1,kk,nx,ny); ind1 = PNTIND(ii,jj,kk,nx,ny);
		} else if (jj == ny-1) {
			ind  = PNTIND(ii,jj,kk,nx,ny);   ind1 = PNTIND(ii,jj-1,kk,nx,ny);
		} else {
			ind  = PNTIND(ii,jj+1,kk,nx,ny); ind1 = PNTIND(ii,jj-1,kk,nx,ny);
			e2v *= 2.;
		}
		for (int ii = 0; ii < f.get_m(); ++ii)
			deri[1 + 3*ii] = (f[ind][ii] - f[ind1][ii])/e2v; // dqdy		

		/* DERIVATIVES WITH RESPECT TO Z */
		if (kk == 0) {
			ind  = PNTIND(ii,jj,kk+1,nx,ny); ind1 = PNTIND(ii,jj,kk,nx,ny);
		} else if (kk == nz-1) {
			ind  = PNTIND(ii,jj,kk,nx,ny);   ind1 = PNTIND(ii,jj,kk-1,nx,ny);
		} else {
			ind  = PNTIND(ii,jj,kk+1,nx,ny); ind1 = PNTIND(ii,jj,kk-1,nx,ny);
			e3w *= 2.;
		}
		for (int ii = 0; ii < f.get_m(); ++ii)
			deri[2 + 3*ii] = (f[ind][ii] - f[ind1][ii])/e3w; // dqdy
	}
	template<class T>
	Field<T> gradOGS2(int nx, int ny, int nz, Field<T> &f, 
		Field<T> &e1, Field<T> &e2, Field<T> &e3, Field<T> &div, Field<T> &curl, Field<T> &Q) {
		// Create output array
		Field<T> grad_UVW(f.get_n(),3*f.get_m(),0.), grad_T(f.get_n(),3*f.get_m());
		// Loop the components
		for (int kk = 0; kk < nz; ++kk) {
			for (int jj = 0; jj < ny; ++jj) { 
				for (int ii = 0; ii < nx; ++ii) {
					// Point id
					int ind = PNTIND(ii,jj,kk,nx,ny);
					int indx = PNTIND(ii-1,jj,kk,nx,ny);
					int indy = PNTIND(ii,jj-1,kk,nx,ny);
					int indz = PNTIND(ii,jj,kk-1,nx,ny);
					// Compute the gradient on UVW
					gradOGS2_ijk(ii,jj,kk,nx,ny,nz,f,e1[ind][1],e2[ind][2],e3[ind][3],grad_UVW[ind]);
					// Project from the UVW grid to T grid
					for (int gg = 0; gg < 3; ++gg) {
						T g1_UVW[3] = {grad_UVW[indx][0 + 3*gg],grad_UVW[indy][1 + 3*gg],grad_UVW[indz][2 + 3*gg]};
						T g2_UVW[3] = {grad_UVW[ind][0 + 3*gg], grad_UVW[ind][1 + 3*gg], grad_UVW[ind][2 + 3*gg]};
						UVW2T_ijk(ii,jj,kk,nx,ny,g1_UVW,g2_UVW,e1,e2,e3,grad_T[ind] + 3*gg);
					}
					// Computation of the divergence
					if (!div.isempty()) {
						div[ind][0] = grad_T[ind][0] + grad_T[ind][4] + grad_T[ind][8];
					}
					// Computation of the curl
					if (!curl.isempty()) {
						curl[ind][0] = grad_T[ind][7] - grad_T[ind][5];
						curl[ind][1] = grad_T[ind][2] - grad_T[ind][6];
						curl[ind][2] = grad_T[ind][3] - grad_T[ind][1];
					}
					// Computation of the Q-criterion
					if (!Q.isempty()) {
						Q[ind][0] = -0.5*(grad_T[ind][0]*grad_T[ind][0] + grad_T[ind][4]*grad_T[ind][4] + grad_T[ind][8]*grad_T[ind][8])
						            -grad_T[ind][1]*grad_T[ind][3]-grad_T[ind][2]*grad_T[ind][6]-grad_T[ind][5]*grad_T[ind][7];
					}
				}
			}
		}
		return grad_T;
	}

	/* GRADOGS4

		Fourth order, face centered, approximation of the gradient for
		a field using OGSTM-BFM approach. This gradient lives in the staggered 
		mesh, and it is projected back to the cell centered mesh. 

		The input is assumed to be on the cell centered mesh.
		This gradient method is experimental.

	*/
	template<class T>
	void gradOGS4_ijk(int ii, int jj, int kk, int nx, int ny, int nz,
		Field<T> &f, T e1u, T e2v, T e3w, T *deri) {

		int ind, ind1, ind2, ind3;

		/* DERIVATIVES WITH RESPECT TO X */
		if (ii <= 1) {
			ind  = PNTIND(ii+2,jj,kk,nx,ny);
			ind1 = PNTIND(ii+1,jj,kk,nx,ny);
			ind2 = PNTIND(ii,jj,kk,nx,ny);

			for (int ii = 0; ii < f.get_m(); ++ii)
				deri[0 + 3*ii] = (-f[ind][ii] + 4.*f[ind1][ii] - 3.*f[ind2][ii])/2./e1u;
		} else if (ii >= nx-2) {
			ind  = PNTIND(ii,jj,kk,nx,ny);
			ind1 = PNTIND(ii-1,jj,kk,nx,ny);
			ind2 = PNTIND(ii-2,jj,kk,nx,ny);

			for (int ii = 0; ii < f.get_m(); ++ii)
				deri[0 + 3*ii] = (3.*f[ind][ii] - 4.*f[ind1][ii] + f[ind2][ii])/2./e1u;
		} else {
			ind  = PNTIND(ii+2,jj,kk,nx,ny);
			ind1 = PNTIND(ii+1,jj,kk,nx,ny);
			ind2 = PNTIND(ii-1,jj,kk,nx,ny);
			ind3 = PNTIND(ii-2,jj,kk,nx,ny);

			for (int ii = 0; ii < f.get_m(); ++ii)
				deri[0 + 3*ii] = (-f[ind][ii] + 8.*f[ind1][ii] - 8.*f[ind2][ii] + f[ind3][ii])/12./e1u;
		}
					
		/* DERIVATIVES WITH RESPECT TO Y */
		if (jj <= 1) {
			ind  = PNTIND(ii,jj+2,kk,nx,ny);
			ind1 = PNTIND(ii,jj+1,kk,nx,ny);
			ind2 = PNTIND(ii,jj,kk,nx,ny);

			for (int ii = 0; ii < f.get_m(); ++ii)
				deri[1 + 3*ii] = (-f[ind][ii] + 4.*f[ind1][ii] - 3.*f[ind2][ii])/2./e2v;
		} else if (jj >= ny-2) {
			ind  = PNTIND(ii,jj,kk,nx,ny);
			ind1 = PNTIND(ii,jj-1,kk,nx,ny);
			ind2 = PNTIND(ii,jj-2,kk,nx,ny);

			for (int ii = 0; ii < f.get_m(); ++ii)
				deri[1 + 3*ii] = (3.*f[ind][ii] - 4.*f[ind1][ii] + f[ind2][ii])/2./e2v;
		} else {
			ind  = PNTIND(ii,jj+2,kk,nx,ny);
			ind1 = PNTIND(ii,jj+1,kk,nx,ny);
			ind2 = PNTIND(ii,jj-1,kk,nx,ny);
			ind3 = PNTIND(ii,jj-2,kk,nx,ny);

			for (int ii = 0; ii < f.get_m(); ++ii)
				deri[1 + 3*ii] = (-f[ind][ii] + 8.*f[ind1][ii] - 8.*f[ind2][ii] + f[ind3][ii])/12./e2v;
		}
	
		/* DERIVATIVES WITH RESPECT TO Z */
		if (kk <= 1) {
			ind  = PNTIND(ii,jj,kk+2,nx,ny);
			ind1 = PNTIND(ii,jj,kk+1,nx,ny);
			ind2 = PNTIND(ii,jj,kk,nx,ny);

			for (int ii = 0; ii < f.get_m(); ++ii)
				deri[2 + 3*ii] = (-f[ind][ii] + 4.*f[ind1][ii] - 3.*f[ind2][ii])/2./e3w;
		} else if (kk >= nz-2) {
			ind  = PNTIND(ii,jj,kk,nx,ny);
			ind1 = PNTIND(ii,jj,kk-1,nx,ny);
			ind2 = PNTIND(ii,jj,kk-2,nx,ny);

			for (int ii = 0; ii < f.get_m(); ++ii)
				deri[2 + 3*ii] = (3.*f[ind][ii] - 4.*f[ind1][ii] + f[ind2][ii])/2./e3w;
		} else {
			ind  = PNTIND(ii,jj,kk+2,nx,ny);
			ind1 = PNTIND(ii,jj,kk+1,nx,ny);
			ind2 = PNTIND(ii,jj,kk-1,nx,ny);
			ind3 = PNTIND(ii,jj,kk-2,nx,ny);

			for (int ii = 0; ii < f.get_m(); ++ii)
				deri[2 + 3*ii] = (-f[ind][ii] + 8.*f[ind1][ii] - 8.*f[ind2][ii] + f[ind3][ii])/12./e3w;
		}
	}
	template<class T>
	Field<T> gradOGS4(int nx, int ny, int nz, Field<T> &f, 
		Field<T> &e1, Field<T> &e2, Field<T> &e3, Field<T> &div, Field<T> &curl, Field<T> &Q) {
		// Create output array
		Field<T> grad_UVW(f.get_n(),3*f.get_m(),0.), grad_T(f.get_n(),3*f.get_m());
		// Loop the components
		for (int kk = 0; kk < nz; ++kk) {
			for (int jj = 0; jj < ny; ++jj) { 
				for (int ii = 0; ii < nx; ++ii) {
					// Point id
					int ind = PNTIND(ii,jj,kk,nx,ny);
					int indx = PNTIND(ii-1,jj,kk,nx,ny);
					int indy = PNTIND(ii,jj-1,kk,nx,ny);
					int indz = PNTIND(ii,jj,kk-1,nx,ny);
					// Compute the gradient on UVW
					gradOGS4_ijk(ii,jj,kk,nx,ny,nz,f,e1[ind][1],e2[ind][2],e3[ind][3],grad_UVW[ind]);
					// Project from the UVW grid to T grid
					for (int gg = 0; gg < 3; ++gg) {
						T g1_UVW[3] = {grad_UVW[indx][0 + 3*gg],grad_UVW[indy][1 + 3*gg],grad_UVW[indz][2 + 3*gg]};
						T g2_UVW[3] = {grad_UVW[ind][0 + 3*gg], grad_UVW[ind][1 + 3*gg], grad_UVW[ind][2 + 3*gg]};
						UVW2T_ijk(ii,jj,kk,nx,ny,g1_UVW,g2_UVW,e1,e2,e3,grad_T[ind] + 3*gg);
					}
					// Computation of the divergence
					if (!div.isempty()) {
						div[ind][0] = grad_T[ind][0] + grad_T[ind][4] + grad_T[ind][8];
					}
					// Computation of the curl
					if (!curl.isempty()) {
						curl[ind][0] = grad_T[ind][7] - grad_T[ind][5];
						curl[ind][1] = grad_T[ind][2] - grad_T[ind][6];
						curl[ind][2] = grad_T[ind][3] - grad_T[ind][1];
					}
					// Computation of the Q-criterion
					if (!Q.isempty()) {
						Q[ind][0] = -0.5*(grad_T[ind][0]*grad_T[ind][0] + grad_T[ind][4]*grad_T[ind][4] + grad_T[ind][8]*grad_T[ind][8])
						            -grad_T[ind][1]*grad_T[ind][3]-grad_T[ind][2]*grad_T[ind][6]-grad_T[ind][5]*grad_T[ind][7];
					}
				}
			}
		}
		return grad_T;
	}

	/* COUNTDEPTHLEVELS

		Counts the number of depth levels (unique values in Z direction) and
		returns the values and the field connectivity.

		If uniquevals is not empty, it will return the depth levels at the
		desired uniquevals.

	*/
	field::Field<int> countDepthLevels(v3::V3v &xyz, std::vector<double> &uniquevals, double epsi);
}

#endif