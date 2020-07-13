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
#include <type_traits>

// Intel Math Kernel Library containing
// optimized BLAS functions
#ifdef MKL
  #include "mkl.h"
#endif

// LAPACK libraries containing
// basic algebra operations
#ifdef USE_LAPACK
	#include "lapacke.h"

	#define LAPACK_MAJOR LAPACK_COL_MAJOR
#endif

#define BLK_LIM 5000

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
			inline matrixMN(const int i)                             { size(i);         fill((T)(0.)); }
			inline matrixMN(const int i, const T val)                { size(i);         fill(val); }
			inline matrixMN(const int i, const T val[])              { size(i);         fill(val); }
			inline matrixMN(const int i, const int j)                { size(i,j);       fill((T)(0.)); }
			inline matrixMN(const int i, const int j, const T val)   { size(i,j);       fill(val); }
			inline matrixMN(const int i, const int j, const T val[]) { size(i,j);       fill(val); }
			inline matrixMN(const matrixMN<T> &mn)                   { size(mn.m,mn.n); fill(mn.mat); }
			inline ~matrixMN()                                       { clear(); }

			// Functions
			inline int  get_m()                                      { return m; }
			inline int  get_n()                                      { return n; }
			inline int  get_bsz()                                    { return bsz; }
			inline void set_bsz(const int i)                         { bsz = i; }
			inline void size(const int i)                            { m = i; n = i; sz = i*i; bsz = -1.; alloc = true; mat = new T[sz]; }
			inline void size(const int i, const int j)               { m = i; n = j; sz = i*j; bsz = -1.; alloc = true; mat = new T[sz]; }
			inline void fill(const T val)                            { if (alloc) std::fill(mat,mat+sz,val); }
			inline void fill(const T val[])                          { if (alloc) std::memcpy(mat, val, sz*sizeof(T)); }
			inline void clear()                                      { n = 0; m = 0; sz = 0; if (alloc) { delete [] mat; } alloc = false; }
			inline T*   data()                                       { return mat; }
			inline void ij(const int i, const int j, const T val)    { mat[n*i + j] = val;  }
			inline T    ij(const int i, const int j) const           { return mat[n*i + j]; }
			inline void print(const char *style); 

			template <class U>
			inline matrixMN<U> convert()                             { matrixMN<U> A(n,m); std::copy(mat,mat+sz,A.data()); return A; }

			// Matrix operations
			inline T            norm2(); // Frobenius norm squared
			inline matrixMN<T>  t();     // Transpose
			inline T            tra();   // Trace
			inline T*           diag();  // Diagonal
			inline void         iden();  // Identity
			inline T            det();   // Determinant
//			inline matrixMN<T>  inv();   // Inverse     (to do)

			// Operations allowed by the LAPACK library 
			#ifdef USE_LAPACK
			// Eigenvalues and eigenvectors
			inline matrixMN<T>      eigen(T wr[], T wi[]);
			inline matrixMN<double> eigen(double wr[], double wi[], matrixMN<double> &V);
			inline matrixMN<float>  eigen(float wr[], float wi[], matrixMN<float> &V);
			// Schur factorization
			inline matrixMN<T> schur(char sort, LAPACK_D_SELECT2 select);
			inline matrixMN<T> schur(char sort, LAPACK_D_SELECT2 select, matrixMN<T> &Z);
			inline matrixMN<double> schur(char sort, LAPACK_D_SELECT2 select, matrixMN<double> &Z, double wr[], double wi[]);
			inline matrixMN<float>  schur(char sort, LAPACK_D_SELECT2 select, matrixMN<float> &Z, float wr[], float wi[]);
			#endif

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
			size_t sz;
			int m, n, bsz;
			T *mat; // Matrix value stored as an array of MxN

			inline matrixMN<T>  t_n();
			inline matrixMN<T>  t_f();
	};

	// Operators
	template<class T>
	inline matrixMN<T>  matrixMN<T>::operator^(const matrixMN<T> &mn) {
		matrixMN<T> out(this->m,this->n);

		#ifdef MKL
			cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,this->m,mn.n,mn.m,1
				this->mat,this->n,mn.mat,1,out.mat,out.n);
		#else
			if (this->m == 2 && this->n == 2) {
				out.ij(0,0, this->ij(0,0)*mn.ij(0,0) + this->ij(0,1)*mn.ij(1,0) );
				out.ij(0,1, this->ij(0,0)*mn.ij(0,1) + this->ij(0,1)*mn.ij(1,1) );
				out.ij(1,0, this->ij(1,0)*mn.ij(0,0) + this->ij(1,1)*mn.ij(1,0) );
				out.ij(1,1, this->ij(1,0)*mn.ij(0,1) + this->ij(1,1)*mn.ij(1,1) );
			} else if (this->m == 3 && this->n == 3) {
				out.ij(0,0, this->ij(0,0)*mn.ij(0,0) + this->ij(0,1)*mn.ij(1,0) + this->ij(0,2)*mn.ij(2,0) );
				out.ij(0,1, this->ij(0,0)*mn.ij(0,1) + this->ij(0,1)*mn.ij(1,1) + this->ij(0,2)*mn.ij(2,1) );
				out.ij(0,2, this->ij(0,0)*mn.ij(0,2) + this->ij(0,1)*mn.ij(1,2) + this->ij(0,2)*mn.ij(2,2) );
				out.ij(1,0, this->ij(1,0)*mn.ij(0,0) + this->ij(1,1)*mn.ij(1,0) + this->ij(1,2)*mn.ij(2,0) );
				out.ij(1,1, this->ij(1,0)*mn.ij(0,1) + this->ij(1,1)*mn.ij(1,1) + this->ij(1,2)*mn.ij(2,1) );
				out.ij(1,2, this->ij(1,0)*mn.ij(0,2) + this->ij(1,1)*mn.ij(1,2) + this->ij(1,2)*mn.ij(2,2) );
				out.ij(2,0, this->ij(2,0)*mn.ij(0,0) + this->ij(2,1)*mn.ij(1,0) + this->ij(2,2)*mn.ij(2,0) );
				out.ij(2,1, this->ij(2,0)*mn.ij(0,1) + this->ij(2,1)*mn.ij(1,1) + this->ij(2,2)*mn.ij(2,1) );
				out.ij(2,2, this->ij(2,0)*mn.ij(0,2) + this->ij(2,1)*mn.ij(1,2) + this->ij(2,2)*mn.ij(2,2) );
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
				for (int j=0; j<this->n; ++j)
					std::printf(style,this->ij(i,j));
				std::printf("|\n");
			}
		} else {
			std::printf("Matrix not allocated!\n");
		}
	}
	template<class T>
	inline T matrixMN<T>::norm2() { // Frobenius norm of a matrix squared.
		T norm2 = 0.;
		for (int i=0; i<this->m; ++i) {
			for (int j=0; j<this->n; ++j)
				norm2 += this->ij(i,j)*this->ij(i,j);
		}
		return norm2;
	}
	template<class T>
	inline matrixMN<T> matrixMN<T>::t_n() { // Naive transpose
		matrixMN<T> out(this->m,this->n); T swp;
		for (int i=0; i<this->m; ++i) {
			for (int j=0; j<i+1; ++j) {
				swp = this->ij(i,j);
				out.ij(i,j,this->ij(j,i));
				out.ij(j,i, swp);
			}
		}
		return out;
	}
	template<class T>
	inline matrixMN<T> matrixMN<T>::t_f() { // Fast transpose (M = N)
		matrixMN<T> out(this->m);
		// Loop by blocks
		for (int ib=0; ib<this->m/this->bsz; ++ib) {
			for(int jb=0; jb<this->n/this->bsz; ++jb) {
				// Loop matrix
				for(int i=ib*this->bsz; i<(ib+1)*this->bsz; ++i) {
					for(int j=jb*this->bsz; j<(jb+1)*this->bsz; ++j) {
						out.ij(j,i, this->ij(i,j)); 
					}
				}// Loop matrix
			}
		}

		return out;
	}	
	template<class T>
	inline matrixMN<T> matrixMN<T>::t() {
		return (this->m == this->n && this-> bsz > 0 && this->m > BLK_LIM) ? this->t_f() : this->t_n();
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

	template<class T>
	inline T matrixMN<T>::det() {
		// Matrix must be n-by-n
		if (this->n != this->m) return 0.;

		T out = 0., sign = 1.;

		if (this->n == 1) {
			out = this->ij(0,0);
		} else if (this->n == 2) {
			out = this->ij(0,0)*this->ij(1,1) - this->ij(0,1)*this->ij(1,0);
		} else {
			matrixMN<T> aux(this->n-1,T(0.));
			for (int p=0; p<this->n; ++p) {
				int subi = 0;
				for (int i=1; i<this->n; ++i) {
					int subj = 0;
					for (int j=0; j<this->n; ++j) {
						if (j == p) continue;
						aux[subi][subj] = this->ij(i,j);
						subj++;
					}
					subi++;
				}
				out += sign*this->ij(0,p)*aux.det();
				sign = -sign;
			}
		}

		return out;
	}

	template<class T>
	inline matrixMN<T> matrixMN<T>::eigen(T wr[], T wi[]) {
		// Matrix to store the eigenvectors
		matrixMN<T> V(this->n);
		// Compute the eigenvalues
		this->eigen(wr,wi,V);
		// Return the eigenvectors as output
		return V;
	}
	template<>
	inline matrixMN<double> matrixMN<double>::eigen(double wr[], double wi[], matrixMN<double> &V) {
		if (this->n != this->m) return (*this); // The matrix is not n-by-n
		// Eigen algorithm overwrites the input matrix, we will copy and return it
		// as the eigenvectors
		matrixMN<double> A(*this);
		double *dummy = new double[this->n];
		// Call lapacke routine
		#ifdef USE_LAPACK
		LAPACKE_dgeev(LAPACK_MAJOR,'N','V',A.get_n(), A.data(), A.get_m(), 
			&wr[0], &wi[0],dummy,this->n,V.data(),this->n);
		#endif
		// Deallocate memory
		delete [] dummy;
		return A; // Returns the factorized matrix
	}
	template<>
	inline matrixMN<float> matrixMN<float>::eigen(float wr[], float wi[], matrixMN<float> &V) {
		if (this->n != this->m) return (*this); // The matrix is not n-by-n
		// Eigen algorithm overwrites the input matrix, we will copy and return it
		// as the eigenvectors
		matrixMN<float> A(*this);
		float *dummy = new float[this->n];
		// Call lapacke routine
		#ifdef USE_LAPACK
		LAPACKE_sgeev(LAPACK_MAJOR,'N','V',A.get_n(), A.data(), A.get_m(), 
			&wr[0], &wi[0],dummy,this->n,V.data(),this->n);
		#endif
		// Deallocate memory
		delete [] dummy;
		return A; // Returns the factorized matrix
	}
	template<class T>
	inline matrixMN<T> matrixMN<T>::schur(char sort, LAPACK_D_SELECT2 select) {
		/*
			The routine computes for an n-by-n real/complex nonsymmetric matrix A, the eigenvalues, 
			the real Schur form T, and, optionally, the matrix of Schur vectors Z. 
			This gives the Schur factorization A = Z*T*ZH.

			Optionally, it also orders the eigenvalues on the diagonal of the real-Schur/Schur form 
			so that selected eigenvalues are at the top left. 
			The leading columns of Z then form an orthonormal basis for the invariant subspace 
			corresponding to the selected eigenvalues.
		*/
		if (this->n != this->m) return (*this); // The matrix is not n-by-n
		// Create the matrix Z
		matrixMN<T> Z(this->n,0.);
		// Run overloaded schur algorithm
		return this->schur(sort,select,Z);
	}

	template<class T>
	inline matrixMN<T> matrixMN<T>::schur(char sort, LAPACK_D_SELECT2 select, matrixMN<T> &Z) {
		/*
			The routine computes for an n-by-n real/complex nonsymmetric matrix A, the eigenvalues, 
			the real Schur form T, and, optionally, the matrix of Schur vectors Z. 
			This gives the Schur factorization A = Z*T*ZH.

			Optionally, it also orders the eigenvalues on the diagonal of the real-Schur/Schur form 
			so that selected eigenvalues are at the top left. 
			The leading columns of Z then form an orthonormal basis for the invariant subspace 
			corresponding to the selected eigenvalues.
		*/
		if (this->n != this->m) return (*this); // The matrix is not n-by-n
		// Preallocate the eigenvalue vectors
		T *wr = new T[this->n], *wi = new T[this->n];
		// Call overloaded schur to obtain the real Schur form
		matrixMN<T> A = this->schur(sort,select,Z,wr,wi);
		// Free memory
		delete [] wr; delete [] wi;
		// Return
		return A;
	}

	template<>
	inline matrixMN<double> matrixMN<double>::schur(char sort, LAPACK_D_SELECT2 select, 
		matrixMN<double> &Z, double wr[], double wi[]) {
		/*
			The routine computes for an n-by-n real/complex nonsymmetric matrix A, the eigenvalues, 
			the real Schur form T, and, optionally, the matrix of Schur vectors Z. 
			This gives the Schur factorization A = Z*T*ZH.

			Optionally, it also orders the eigenvalues on the diagonal of the real-Schur/Schur form 
			so that selected eigenvalintues are at the top left. 
			The leading columns of Z then form an orthonormal basis for the invariant subspace 
			corresponding to the selected eigenvalues.
		*/
		if (this->n != this->m) return (*this); // The matrix is not n-by-n
		// Schur algorithm overwrites the input matrix, we will copy and return it
		// as the schur transformation
		int sdim = 0; matrixMN<double> A(*this);
		// wr, wi and Z should already be allocated at this point
		// Call lapacke routine
		#ifdef USE_LAPACK
		LAPACKE_dgees(LAPACK_MAJOR, 'V', sort, select, A.get_n(), A.data(), 
			A.get_n(), &sdim, &wr[0], &wi[0], Z.data(), A.get_n());
		#endif
		// Return A as the real Schur form
		return A;		
	}

	template<>
	inline matrixMN<float> matrixMN<float>::schur(char sort, LAPACK_D_SELECT2 select, 
		matrixMN<float> &Z, float wr[], float wi[]) {
		/*
			The routine computes for an n-by-n real/complex nonsymmetric matrix A, the eigenvalues, 
			the real Schur form T, and, optionally, the matrix of Schur vectors Z. 
			This gives the Schur factorization A = Z*T*ZH.

			Optionally, it also orders the eigenvalues on the diagonal of the real-Schur/Schur form 
			so that selected eigenvalues are at the top left. 
			The leading columns of Z then form an orthonormal basis for the invariant subspace 
			corresponding to the selected eigenvalues.
		*/
		if (this->n != this->m) return (*this); // The matrix is not n-by-n
		// Schur algorithm overwrites the input matrix, we will copy and return it
		// as the schur transformation
		int sdim = 0; matrixMN<float> A(*this);
		#ifdef USE_LAPACK
		// wr, wi and Z should already be allocated at this point
		// Lapack works in doubles
		matrixMN<double> A_aux = A.convert<double>(), Z_aux = Z.convert<double>();
		double *wr_aux = new double[this->n], *wi_aux = new double[this->n];
		// Call lapacke routine
		LAPACKE_dgees(LAPACK_MAJOR, 'V', sort, select, A_aux.get_n(), A_aux.data(), 
			A_aux.get_n(), &sdim, &wr_aux[0], &wi_aux[0], Z_aux.data(), A_aux.get_n());
		// Copy outputs
		std::copy(wr_aux,wr_aux+this->n,wr); std::copy(wi_aux,wi_aux+this->n,wi);
		A = A_aux.convert<float>(); Z = Z_aux.convert<float>();
		delete [] wr_aux; delete [] wi_aux;
		#endif
		// Return A as the real Schur form
		return A;		
	}
	#endif
}
#endif