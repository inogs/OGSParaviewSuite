/*
	MATRIX(M,N)

	Defines a 2 dimensional matrix of dimensions M x N. 
	Allows for operations between matrices.

	Let A i B be matrices and a a numeric value...
	
	Usage:
		> A = B sets the contents of B to A if size(A) == size(B)
		> A + a sums a to all positions of A
		> A - a substracts a to all positions of A
		> A * a product a to all positions of A
		> A / a division a to all positions of A
		> A + B sums A and B if size(A) == size(B)
		> A - B substracts A and B if size(A) == size(B)
		> A * B point to point product of A and B
		> A / B point to point division of A and B
		> A ^ B matrix dot product
		> A == B is testing size(A) == size(B), contents may be different!!
		> A != B is testing size(A) != size(B), contents may be different!!

    Functions:
     	> size(i) or size(i,j) defines the size of the matrix
     	> ij(i,j) or ij(i,j,val) lets grab a value using i,j indexes (from 1 to end)
     	  and also set a new value (val) at the i,j position
     	> print(syle) prints the matrix in a nice printf style

	Arnau Miro (UPC-ESEIAAT) (c) 2018
*/

#ifndef matrixMN_h
#define matrixMN_h

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <algorithm>

// Intel Math Kernel Library containing
// optimized BLAS functions
#ifdef MKL
  #include "mkl.h"
#endif

namespace matMN
{
	template<class T>
	class matrixMN {

		friend matrixMN<T> operator+(T val, const matrixMN<T> &mn) { matrixMN<T> out(mn.m,mn.n); for (int i=0; i<mn.sz; ++i) {out.mat[i] = val + mn.mat[i];} return out; }
		friend matrixMN<T> operator-(T val, const matrixMN<T> &mn) { matrixMN<T> out(mn.m,mn.n); for (int i=0; i<mn.sz; ++i) {out.mat[i] = val - mn.mat[i];} return out; }
		friend matrixMN<T> operator*(T val, const matrixMN<T> &mn) { matrixMN<T> out(mn.m,mn.n); for (int i=0; i<mn.sz; ++i) {out.mat[i] = val * mn.mat[i];} return out; }
		friend matrixMN<T> operator/(T val, const matrixMN<T> &mn) { matrixMN<T> out(mn.m,mn.n); for (int i=0; i<mn.sz; ++i) {out.mat[i] = val / mn.mat[i];} return out; }
	
		public:
			// Constructors and destructor
			inline matrixMN()                                        { m = 0; n = 0; sz = 0; }
			inline matrixMN(const int i)                             { size(i);         fill(0.); }
			inline matrixMN(const int i, const T val)                { size(i);         fill(val); }
			inline matrixMN(const int i, const T val[])              { size(i);         fill(val); }
			inline matrixMN(const int i, const int j)                { size(i,j);       fill(0.); }
			inline matrixMN(const int i, const int j, const T val)   { size(i,j);       fill(val); }
			inline matrixMN(const int i, const int j, const T val[]) { size(i,j);       fill(val); }
			inline matrixMN(const matrixMN<T> &mn)                   { size(mn.m,mn.n); fill(mn.mat); }
			inline ~matrixMN()                                       { clear(); }

			// Functions
			inline int  get_m()                                      { return m; }
			inline int  get_n()                                      { return n; }
			inline void size(const int i)                            { m = i; n = i; sz = i*i; alloc = true; mat = new T[sz]; }
			inline void size(const int i, const int j)               { m = i; n = j; sz = i*j; alloc = true; mat = new T[sz]; }
			inline void fill(const T val)                            { if (alloc) std::fill(mat,mat+sz,val); }
			inline void fill(const T val[])                          { if (alloc) std::memcpy(mat, val, sz*sizeof(T)); }
			inline void clear()                                      { n = 0; m = 0; sz = 0; if (alloc) { delete [] mat; } alloc = false; }
			inline T*   data()                                       { return mat; }
			inline void ij(const int i, const int j, const T val)    { mat[n*i + j] = val;  }
			inline T    ij(const int i, const int j) const           { return mat[n*i + j]; }
			inline void print(const char *style); 

			// Matrix operations
			inline T            norm2(); // Frobenius norm squared
			inline matrixMN<T>  t();     // Transpose
			inline T            tra();   // Trace
			inline T*           diag();  // Diagonal
			inline void         iden();  // Identity
//			inline matrixMN<T>  inv();   // Inverse (to do)
//			inline T            det();   // Determinant (to do)

			// Operators
			inline matrixMN<T> &operator=(const matrixMN<T> &mn)        { clear(); size(mn.m,mn.n); fill(mn.mat); return (*this);}                                                            // Equality of matrix sizes
			inline matrixMN<T>  operator+(const T val) const            { matrixMN<T> out(m,n); for (int i=0; i<sz; ++i) {out.mat[i] = mat[i] + val;} return out; }                           // Sum matrix and number at all points
			inline matrixMN<T>  operator-(const T val) const            { matrixMN<T> out(m,n); for (int i=0; i<sz; ++i) {out.mat[i] = mat[i] - val;} return out; }                           // Substract matrix and number at all points
			inline matrixMN<T>  operator*(const T val) const            { matrixMN<T> out(m,n); for (int i=0; i<sz; ++i) {out.mat[i] = mat[i] * val;} return out; }                           // Product matrix and number at all points
			inline matrixMN<T>  operator/(const T val) const            { matrixMN<T> out(m,n); for (int i=0; i<sz; ++i) {out.mat[i] = mat[i] / val;} return out; }                           // Division matrix and number at all points
			inline matrixMN<T>  operator+(const matrixMN<T> &mn)        { matrixMN<T> out(m,n); if ((*this) == mn) { for(int i=0; i<sz; ++i){out.mat[i] = mat[i] + mn.mat[i];}} return out; } // Sum of two matrices
			inline matrixMN<T>  operator-(const matrixMN<T> &mn)        { matrixMN<T> out(m,n); if ((*this) == mn) { for(int i=0; i<sz; ++i){out.mat[i] = mat[i] - mn.mat[i];}} return out; } // Substraction of two matrices
			inline matrixMN<T>  operator*(const matrixMN<T> &mn)        { matrixMN<T> out(m,n); if ((*this) == mn) { for(int i=0; i<sz; ++i){out.mat[i] = mat[i] * mn.mat[i];}} return out; } // Point to point product
			inline matrixMN<T>  operator/(const matrixMN<T> &mn)        { matrixMN<T> out(m,n); if ((*this) == mn) { for(int i=0; i<sz; ++i){out.mat[i] = mat[i] / mn.mat[i];}} return out; } // Point to point division
			inline matrixMN<T>  operator^(const matrixMN<T> &mn);                                                                                                                             // Matrix dot product (to do)
			inline void         operator+=(const T val)                 { for (int i=0; i<sz; ++i) mat[i] += val; }                                                                           // Sum matrix and number at all points
			inline void         operator-=(const T val)                 { for (int i=0; i<sz; ++i) mat[i] -= val; }                                                                           // Substract matrix and number at all points
			inline void         operator*=(const T val)                 { for (int i=0; i<sz; ++i) mat[i] *= val; }                                                                           // Product matrix and number at all points
			inline void         operator/=(const T val)                 { for (int i=0; i<sz; ++i) mat[i] /= val; }                                                                           // Division matrix and number at all points
			inline void         operator+=(const matrixMN<T> &mn)       { if ((*this) == mn) {for(int i=0; i<sz; ++i) {mat[i] += mn.mat[i];}} }                                               // Sum of two matrices
			inline void         operator-=(const matrixMN<T> &mn)       { if ((*this) == mn) {for(int i=0; i<sz; ++i) {mat[i] -= mn.mat[i];}} }                                               // Substraction of two matrices
			inline void         operator*=(const matrixMN<T> &mn)       { if ((*this) == mn) {for(int i=0; i<sz; ++i) {mat[i] *= mn.mat[i];}} }                                               // Point to point product
			inline void         operator/=(const matrixMN<T> &mn)       { if ((*this) == mn) {for(int i=0; i<sz; ++i) {mat[i] /= mn.mat[i];}} }                                               // Point to point division
			inline bool         operator==(const matrixMN<T> &mn) const { return( (m == mn.m && n == mn.n) ? true : false ); }                                                                // This is used to test equality of sizes
			inline T*           operator[](int i)                       { return mat + n*i; }
			
		private:
			bool alloc = false;
			int m, n, sz;
			T *mat; // Matrix value stored as an array of MxN
	};

	// Operators
	template<class T>
	inline matrixMN<T>  matrixMN<T>::operator^(const matrixMN<T> &mn) {
		matrixMN<T> out(this->m,this->n);

		#ifdef MKL
			cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,this->m,mn.n,mn.m,1
				this->mat,this->n,mn.mat,1,out.mat,out.n);
		#else
			if (this->m = 3 && this->n == 3) {
				out.mat[0] = this->mat[0]*mn.mat[0] + this->mat[1]*mn.mat[3] + this->mat[2]*mn.mat[6];
				out.mat[1] = this->mat[0]*mn.mat[1] + this->mat[1]*mn.mat[4] + this->mat[2]*mn.mat[7];
				out.mat[2] = this->mat[0]*mn.mat[2] + this->mat[1]*mn.mat[5] + this->mat[2]*mn.mat[8];
				out.mat[3] = this->mat[3]*mn.mat[0] + this->mat[4]*mn.mat[3] + this->mat[5]*mn.mat[6];
				out.mat[4] = this->mat[3]*mn.mat[1] + this->mat[4]*mn.mat[4] + this->mat[5]*mn.mat[7];
				out.mat[5] = this->mat[3]*mn.mat[2] + this->mat[4]*mn.mat[5] + this->mat[5]*mn.mat[8];
				out.mat[6] = this->mat[6]*mn.mat[0] + this->mat[7]*mn.mat[3] + this->mat[8]*mn.mat[6];
				out.mat[7] = this->mat[6]*mn.mat[1] + this->mat[7]*mn.mat[4] + this->mat[8]*mn.mat[7];
				out.mat[8] = this->mat[6]*mn.mat[2] + this->mat[7]*mn.mat[5] + this->mat[8]*mn.mat[8];
			} else {
				for (int ii=0; ii<out.m; ii++) {
					for (int jj=0; jj<out.n; jj++) {
						for(int kk=0; kk<this->n; kk++)
							out.mat[this->n*ii+jj] += this->ij(ii,kk)*mn.ij(kk,ii);
					}
				}
			}
		#endif

		return out;
	}	

	// Functions
	template<class T>
	inline void matrixMN<T>::print(const char *style) {
		if (alloc) {
			for (int i=0; i<this->m; ++i) {
				std::printf("|");
				for (int j=0; j<this->n; j++)
					std::printf(style,this->ij(i,j));
				std::printf("|\n");
			}
		} else {
			std::printf("Matrix not allocated!\n");
		}
	}
	template<class T>
	inline T matrixMN<T>::norm2() {
		// Frobenius norm of a matrix squared.
		T norm2 = 0.;
		for (int i=0; i<this->m; ++i) {
			for (int j=0; j<this->n; j++)
				norm2 += this->ij(i,j)*this->ij(i,j);
		}
		return norm2;
	}
	template<class T>
	inline matrixMN<T> matrixMN<T>::t() {
		matrixMN<T> out(this->m,this->n);
		for (int i=0; i<this->m; ++i) {
			for (int j=0; j<this->n; j++)
				out.ij(i,j,this->ij(j,i));
		}
		return out;
	}
	template<class T>
	inline T *matrixMN<T>::diag() {
		T *diag = NULL; 
		if (this->m == this->n) {
			diag = new T[this->n];
			for (int i=0; i<this->m; ++i)
				diag[i] = this->ij(i,i);
		}
		return diag;
	}
	template<class T>
	inline T matrixMN<T>::tra() {
		T tra = 0.;
		if (this->m == this->n) {
			for (int i=0; i<this->m; ++i)
				tra += this->ij(i,i);
		}
		return tra;
	}
	template<class T>
	inline void matrixMN<T>::iden() {
		this->fill(0.);
		if (this->m == this->n) {
			for (int i=0; i<this->m; ++i)
				this->ij(i,i,1.);
		}
	}

/* EIGEN carries out the Jacobi eigenvalue iteration.
 *
 * This function computes the eigenvalues and eigenvectors of a
 * real symmetric matrix, using Rutishauser's modfications of the classical
 * Jacobi rotation method with threshold pivoting. 
 *
/*----------------------------------------------------------------------------*/
	template<class T>
	int eigen(matrixMN<T> &mat, T d[], matrixMN<T> &v, int maxiter) {

		// Initalize eigenvectors and eigenvalues
		v.iden();     // Initialize v to the identity matrix

		// Extra variables
		int n_iter = 0, n_rot = 0;
		
		T *bw = mat.diag(); std::memcpy(d,bw,mat.get_n()*sizeof(T));

		T *zw; zw = new T[mat.get_n()]; std::fill(zw,zw+v.get_n(),0.);

		// Loop iterations
		for (n_iter = 0; n_iter < maxiter; n_iter++) {
			// The convergence threshold is based on the size of the elements in
			// the strict upper triangle of the matrix.
			double thresh = 0.;
			for (int ii=0; ii<mat.get_n(); ++ii)
				for (int jj=0; jj<ii; ++jj)
					thresh += mat.ij(ii,jj)*mat.ij(ii,jj);
				thresh = std::sqrt(thresh)/(double)(4.*v.get_n());

			if (thresh == 0.) break;

			// Jacobi algorithm
			for (int p=0; p<mat.get_n(); p++) {
				for (int q=p+1; q<mat.get_n(); q++) {
					double gapq  = 10.0 * std::fabs(mat.ij(q,p));
					double termp = gapq + std::fabs(d[p]);
					double termq = gapq + std::fabs(d[q]);
					//Annihilate tiny offdiagonal elements
					if (4<n_iter && termp == std::fabs(d[p]) && termq == std::fabs(d[q]) )
						mat.ij(p,q,0.);
					// Otherwise, apply a rotation.
					else if ( thresh <= std::fabs(mat.ij(q,p)) ) { 
						double h    = d[q] - d[p];
						double term = std::fabs(h) + gapq;

						double t;
						if (term == std::fabs(h))
							t = mat.ij(q,p)/h;
						else {
							double theta = 0.5*h/mat.ij(q,p);
							t = 1.0/( std::fabs(theta) + std::sqrt(1.0 + theta*theta) );
							t = (theta < 0.) ? -t : t;
						}

						double c   = 1.0/std::sqrt(1.0 + t*t);
						double s   = t*c;
						double tau = s/(1.0 + c);

						h = t * mat.ij(q,p);

						//Accumulate corrections to diagonal elements.
						zw[p] -= h;                 
						zw[q] += h;
						d[p]  -= h;
						d[q]  += h;

						mat.ij(q,p,0.);

						//  Rotate, using information from the upper triangle of A only.
						double g;
						for (int j=0; j<p; ++j) {
							g = mat.ij(p,j);
							h = mat.ij(q,j);
							mat.ij(p,j, g - s*(h + g*tau) );
							mat.ij(q,j, h + s*(g - h*tau) );
						}

						for (int j=p+1; j<q; ++j) {
							g = mat.ij(j,p);
							h = mat.ij(q,j);
							mat.ij(j,p, g - s*(h + g*tau) );
							mat.ij(q,j, h + s*(g - h*tau) );
						}

						for (int j=q+1; j<mat.get_n(); ++j) {
							g = mat.ij(j,p);
							h = mat.ij(j,q);
							mat.ij(j,p, g - s*(h + g*tau) );
							mat.ij(j,q, h + s*(g - h*tau) );
						}

						//Accumulate information in the eigenvector matrix.
						for (int j=0; j<mat.get_n(); ++j) {
							g = v.ij(p,j);
							h = v.ij(q,j);
							v.ij(p,j, g - s*(h + g*tau) );
							v.ij(q,j, h + s*(g - h*tau) );
						}

						++n_rot;
					}
				}
			}

			for (int i=0; i<mat.get_n(); ++i) {
				bw[i] = bw[i] + zw[i];
				d[i]  = bw[i];
				zw[i] = 0.0;
			}
		}

		// Restore lower triangle of input matrix.
		for (int i = 0; i < mat.get_n(); ++i)
			for (int j = 0; j < i; j++ )
				mat.ij(i,j, mat.ij(j,i) );

		//  Ascending sort the eigenvalues and eigenvectors.
		for (int k=0; k<mat.get_n()-1; k++) {
			int m = k;
			for (int l=k+1; l<mat.get_n(); ++l)
				m = (d[l] < d[m]) ? l : m;

			if (m != k) {
				T t = d[m];
				d[m]     = d[k];
				d[k]     = t;

				for (int i=0; i<mat.get_n(); ++i) {
					T w = v.ij(m,i);
					v.ij(m,i,v.ij(k,i));
					v.ij(k,i,w);
				}
			}
		}
		delete [] bw; delete [] zw;
		return n_iter;
	}
}

#endif