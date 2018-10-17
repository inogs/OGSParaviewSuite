/*=========================================================================

  Program:   Utilities
  Module:    matrixOperations.cxx

  Useful matrix operations extracted from Code_Saturne and the web.

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include <cmath>

#define MATH_EPZERO 1.e-12
#define MATH_PI     3.141592653589793

/*----------------------------------------------------------------------------*/
/* Compute the product of a matrix of 3x3 real values by a matrix of 3x3
 * real values.
 *
 * \param[in]     A             matrix of 3x3 real values
 * \param[in]     B             matrix of 3x3 real values
 * \param[out]    C             matrix of 3x3 real values
 */
/*----------------------------------------------------------------------------*/
static inline void math_33_3_product(const double A[9], 
                                     const double B[9], 
                                     double *C) {
  C[0] = A[0]*B[0] + A[1]*B[3] + A[2]*B[6];
  C[1] = A[0]*B[1] + A[1]*B[4] + A[2]*B[7];
  C[2] = A[0]*B[2] + A[1]*B[5] + A[2]*B[8];
  C[3] = A[3]*B[0] + A[4]*B[3] + A[5]*B[6];
  C[4] = A[3]*B[1] + A[4]*B[4] + A[5]*B[7];
  C[5] = A[3]*B[2] + A[4]*B[5] + A[5]*B[8];
  C[6] = A[6]*B[0] + A[7]*B[3] + A[8]*B[6];
  C[7] = A[6]*B[1] + A[7]*B[4] + A[8]*B[7];
  C[8] = A[6]*B[2] + A[7]*B[5] + A[8]*B[8];
}

static inline void math_33_3_product_add(const double A[9], 
                                         const double B[9], 
                                         double *C) {
  C[0] += A[0]*B[0] + A[1]*B[3] + A[2]*B[6];
  C[1] += A[0]*B[1] + A[1]*B[4] + A[2]*B[7];
  C[2] += A[0]*B[2] + A[1]*B[5] + A[2]*B[8];
  C[3] += A[3]*B[0] + A[4]*B[3] + A[5]*B[6];
  C[4] += A[3]*B[1] + A[4]*B[4] + A[5]*B[7];
  C[5] += A[3]*B[2] + A[4]*B[5] + A[5]*B[8];
  C[6] += A[6]*B[0] + A[7]*B[3] + A[8]*B[6];
  C[7] += A[6]*B[1] + A[7]*B[4] + A[8]*B[7];
  C[8] += A[6]*B[2] + A[7]*B[5] + A[8]*B[8];
}

/*----------------------------------------------------------------------------*/
/* Compute the square norm of a vector of 3 real values.
 *
 * \param[in]     v             vector of 3 real values
 *
 * \return square norm of v.
 *
/*----------------------------------------------------------------------------*/
static inline double math_3_square_norm(const double v[3]) { 
	return (v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

/*----------------------------------------------------------------------------*/
/* Compute the determinant of a 3x3 symmetric matrix
 *
 * \param[in]  m    3x3 symmetric matrix
 *
 * \return the determinant
 */
/*----------------------------------------------------------------------------*/
static inline double math_sym_33_determinant(const double m[6]) {
  const double com0 = m[1]*m[2] - m[4]*m[4];
  const double com1 = m[4]*m[5] - m[3]*m[2];
  const double com2 = m[3]*m[4] - m[1]*m[5];

  return m[0]*com0 + m[3]*com1 + m[5]*com2;
}

/*----------------------------------------------------------------------------*/
/* Compute all eigenvalues of a 3x3 symmetric matrix with symmetric storage.
 *
 * Based on: Oliver K. Smith "eigenvalues of a symmetric 3x3 matrix",
 *           Communication of the ACM (April 1961)
 *           (Wikipedia article entitled "Eigenvalue algorithm")
 *
 * \param[in]  m          3x3 symmetric matrix (m11, m22, m33, m12, m23, m13)
 * \param[out] eig_vals   size 3 vector
 *
/*----------------------------------------------------------------------------*/
void math_sym_33_eigen(const double m[6], double eig_vals[3]) {
  double  e, e1, e2, e3;

  double  p1 = math_3_square_norm((const double *)(m+3));
  double  d2 = math_3_square_norm((const double *)m);

  if (p1 > MATH_EPZERO*d2) { /* m is not diagonal */

    double  n[6];
    double  tr = (m[0] + m[1] + m[2]);
    double  tr_third = 1./3. * tr;

    e1 = m[0] - tr_third, e2 = m[1] - tr_third, e3 = m[2] - tr_third;
    double  p2 = e1*e1 + e2*e2 + e3*e3 + 2.*p1;

    double  p = sqrt(p2*1./6.);
    double  ovp = 1./p;

    for (int  i = 0; i < 3; i++) {
      /* Diagonal */
      n[i] = ovp * (m[i] - tr_third);
      /* Extra diagonal */
      n[3 + i] = ovp * m[3 + i];
    }

    /* r should be between -1 and 1 but truncation error and bad conditionning
       can lead to slighty under/over-shoot */
    double  r = 0.5 * math_sym_33_determinant(n);

    double  cos_theta, cos_theta_2pi3;
    if (r <= -1.) {
      cos_theta = 0.5; // theta = pi/3;
      cos_theta_2pi3 = -1.;
    }
    else if (r >= 1.) {
      cos_theta = 1.; // theta = 0.;
      cos_theta_2pi3 = -0.5;
    }
    else {
      cos_theta = cos(1./3.*acos(r));
      cos_theta_2pi3 = cos(1./3.*(acos(r) + 2.*MATH_PI));
    }

    /* eigenvalues computed should satisfy e1 < e2 < e3 */
    e3 = tr_third + 2.*p*cos_theta;
    e1 = tr_third + 2.*p*cos_theta_2pi3;
    e2 = tr - e1 -e3; // since tr(m) = e1 + e2 + e3

  }
  else { // m is diagonal

    e1 = m[0], e2 = m[1], e3 = m[2];

  } /* diagonal or not */

  if (e3 < e2) e = e3, e3 = e2, e2 = e;
  if (e3 < e1) e = e3, e3 = e1, e1 = e2, e2 = e;
  else {
    if (e2 < e1) e = e2, e2 = e1, e1 = e;
  }
  /* Return values */
  eig_vals[0] = e1;
  eig_vals[1] = e2;
  eig_vals[2] = e3;
}

/*----------------------------------------------------------------------------*/
/* MATH_MAT_DIAG_GET_VECTOR gets the value of the diagonal of a matrix.  
 *
 * The matrix is a doubly dimensioned array of double values, stored as a vector
 * in row-major order.
 *
 * Input, int N, the number of rows and columns of the matrix.
 * Input, double A[N*N], the N by N matrix.
 * Output, double V[N], the diagonal entries of the matrix.
 *
/*----------------------------------------------------------------------------*/
inline void math_mat_diag_get_vector(int n, double a[], double v[]) {
  for (int i = 0; i < n; i++) 
    v[i] = a[n*i+i];
}

/*----------------------------------------------------------------------------*/
/* MATH_MAT_IDENTITY sets the square matrix A to the identity.
 *
 * The matrix is a doubly dimensioned array of double values, stored as a vector 
 * in row-major order.
 *
 * Input, int N, the order of A.
 * Output, double A[N*N], the N by N identity matrix.
 *
/*----------------------------------------------------------------------------*/
inline void math_mat_identity (int n, double a[]) {
  int k = 0;
  for (int j = 0; j < n; j++ ) {
    for (int i = 0; i < n; i++ ) {
      a[k] = (i == j) ? 1.0 : 0.0;
      k++;
    }
  }
}

/*----------------------------------------------------------------------------*/
/* MATH_MAT_NORM_FRO returns the Frobenius norm of a matrix.
 *
 * The matrix is a doubly dimensioned array of double values, stored as a vector 
 * in row-major order.
 *
 * The Frobenius norm is defined as
 *
 *    NORM_FRO = sqrt( sum ( 1 <= I <= M ) sum ( 1 <= j <= N ) A(I,J)^2 )
 *
 * The matrix Frobenius norm is not derived from a vector norm, but
 * is compatible with the vector L2 norm, so that:
 *
 *    norm_l2 ( A * x ) <= r8mat_norm_fro ( A ) * r8vec_norm_l2 ( x ).
 *
 * Input, int N, the order of A.
 * Output, double A[N*N], the N by N identity matrix.
 *
/*----------------------------------------------------------------------------*/
inline double math_mat_norm_fro (int m, int n, double a[]) {
  double value = 0.;

  for (int j = 0; j < n; j++) {
    for (int i = 0; i < m; i++ ) {
      value += a[n*i + j]*a[n*i + j];
    }
  }
  value = sqrt(value);

  return value;
}

/*----------------------------------------------------------------------------*/
/* MATH_JACOBI_EIGENVALUE carries out the Jacobi eigenvalue iteration.
 *
 * This function computes the eigenvalues and eigenvectors of a
 * real symmetric matrix, using Rutishauser's modfications of the classical
 * Jacobi rotation method with threshold pivoting. 
 *
 * Input, int N, the order of the matrix.
 * Input, double A[N*N], the matrix, which must be square, real, and symmetric.
 * Input, int IT_MAX, the maximum number of iterations.
 * Output, double V[N*N], the matrix of eigenvectors.
 * Output, double D[N], the eigenvalues, in descending order.
 * Output, int &IT_NUM, the total number of iterations.
 * Output, int &ROT_NUM, the total number of rotations.
 *
/*----------------------------------------------------------------------------*/

void math_jacobi_eigenvalue (int n, double a[], int it_max, double v[], 
  double d[], int &it_num, int &rot_num ) {
  
  double *bw, *zw;

  bw = new double[n];
  zw = new double[n];

  math_mat_identity(n,v);
  math_mat_diag_get_vector(n,a,d);

  for (int i = 0; i < n; i++ ) {
    bw[i] = d[i];
    zw[i] = 0.0;
  }
  it_num = 0; rot_num = 0;

  while ( it_num < it_max ) {
    it_num++;

    //The convergence threshold is based on the size of the elements in
    //the strict upper triangle of the matrix.
    double thresh = 0.0;
    for (int j = 0; j < n; j++ ) {
      for (int i = 0; i < j; i++ ) {
        thresh += a[n*i+j] * a[n*i+j];
      }
    }
    thresh = sqrt(thresh)/(double)(4.*n);

    if (thresh == 0.0) break;

    double gapq, termp, termq;
    for (int p = 0; p < n; p++ ) {
      for (int q = p + 1; q < n; q++ ) {
        gapq  = 10.0 * fabs ( a[n*p+q] );
        termp = gapq + fabs ( d[p] );
        termq = gapq + fabs ( d[q] );
        
        //  Annihilate tiny offdiagonal elements.
        if ( 4 < it_num && termp == fabs ( d[p] ) && termq == fabs ( d[q] ) ) {
          a[n*p+q] = 0.0;
        } else if ( thresh <= fabs ( a[p+q*n] ) ) {// Otherwise, apply a rotation.
          double h    = d[q] - d[p];
          double term = fabs ( h ) + gapq;

          double t, theta;
          if ( term == fabs ( h ) ) {
            t = a[p+q*n] / h;
          } else {
            theta = 0.5 * h / a[n*p+q];
            t     = 1.0 / ( fabs ( theta ) + sqrt ( 1.0 + theta * theta ) );
            if ( theta < 0.0 ) t = - t;
          }
          double c   = 1.0 / sqrt ( 1.0 + t * t );
          double s   = t * c;
          double tau = s / ( 1.0 + c );
          h = t * a[p+q*n];

          // Accumulate corrections to diagonal elements.
          zw[p] = zw[p] - h;                 
          zw[q] = zw[q] + h;
          d[p]  =  d[p] - h;
          d[q]  =  d[q] + h;

          a[n*p+q] = 0.0;

          // Rotate, using information from the upper triangle of A only.
          double g;
          for (int j = 0; j < p; j++ ) {
            g        = a[n*j+p];
            h        = a[n*j+q];
            a[n*j+p] = g - s * ( h + g * tau );
            a[n*j+q] = h + s * ( g - h * tau );
          }

          for (int j = p + 1; j < q; j++ ) {
            g        = a[n*p+j];
            h        = a[n*j+q];
            a[n*p+j] = g - s * ( h + g * tau );
            a[n*j+q] = h + s * ( g - h * tau );
          }

          for (int j = q + 1; j < n; j++ ) {
            g        = a[n*p+j];
            h        = a[n*q+j];
            a[n*p+j] = g - s * ( h + g * tau );
            a[n*q+j] = h + s * ( g - h * tau );
          }

          // Accumulate information in the eigenvector matrix.
          for (int j = 0; j < n; j++ ) {
            g        = v[n*j+p];
            h        = v[n*j+q];
            v[n*j+p] = g - s * ( h + g * tau );
            v[n*j+q] = h + s * ( g - h * tau );
          }
          rot_num++;
        }
      }
    }

    for (int i = 0; i < n; i++ ) {
      bw[i] = bw[i] + zw[i];
      d[i]  = bw[i];
      zw[i] = 0.0;
    }
  }

  // Restore upper triangle of input matrix.
  for (int j = 0; j < n; j++ ) {
    for (int i = 0; i < j; i++ ) {
      a[n*i+j] = a[n*j+i];
    }
  }

  // Ascending sort the eigenvalues and eigenvectors.
  for (int k = 0; k < n - 1; k++ ) {
    int m = k;
    for (int l = k + 1; l < n; l++ ) {
      if ( d[l] < d[m] ) m = l;
    }

    if ( m != k ) {
      double t = d[m];
      d[m] = d[k];
      d[k] = t;
      for (int i = 0; i < n; i++ ) {
        double w = v[n*i+m];
        v[n*i+m] = v[n*i+k];
        v[n*i+k] = w;
      }
    }
  }
  delete [] bw, zw;
}