/* 
    VECTOR 3D
	
	Defines a 3 dimensional vector composed of |x,y,z|.
	Allows operations between vectors.

	Let v1 and v2 be vectors of |x,y,z| and k a random value...

	Construction:
		> V3()      empty construction (sets to zero)
		> V3(x,y,z) build from a series of x,y,z points
		> V3(*x)    build from a C/C++ array
		> V3( V3 )  build from an existing V3 vector
	
	Usage:
		> v1 = v2   sets the contents of v2 to v1
		> k  + v1   adds a V3 vector by k in all positions
		> v1 + k    adds a V3 vector by k in all positions
		> k  - v1   substracts a V3 vector by k in all positions
		> v1 - k    substracts a V3 vector by k in all positions
		> k  * v1   multiplies a V3 vector by k in all positions
		> v1 * k    multiplies a V3 vector by k in all positions
		> v1 / k    divides a V3 vector by k in all positions
		> v1 + v2   adds two V3 vectors
		> v1 - v2   subtracts two V3 vectors
		> v1 * v2   scalar product between two V3 vectors
		> v1 ^ v2   cross product between two V3 vectors
		> v1 += k   adds a V3 vector by k in all positions
		> v1 -= k   substracts a V3 vector by k in all positions
		> v1 *= k   multiplies a V3 vector by k in all positions
		> v1 /= k   divides a V3 vector by k in all positions
		> v1 += v2  adds two V3 vectors
		> v1 -= v2  subtracts two V3 vectors
		> v1 *= v2  point to point product between two V3 vectors
		> v1 == v2  checks equalty using the euclidean norm

	Functions:
		> norm2:    performs the square of the vector norm
		> print:    prints the vector in a nice format

	Arnau Miro (UPC-ESEIAAT) (c) 2018
*/

#ifndef V3_h
#define V3_h

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <algorithm>

namespace v3
{
	class V3 {

		friend inline V3 operator+(double val, const V3 &v);
		friend inline V3 operator-(double val, const V3 &v);
		friend inline V3 operator*(double val, const V3 &v);
		
		public:
			// Constructor
			inline V3();
			inline V3(const double x, const double y, const double z);
			inline V3(const V3 &v);
			inline V3(const double *x);
			// Destructor
			inline ~V3();
	
			// Functions
			inline double norm2() const;
			inline void   print();

			// Operators
			inline V3     &operator=(const V3 &v);

			inline V3      operator+(const double v) const; // Sum of a number to a vector
			inline V3      operator-(const double v) const; // Substraction of a number to a vector
			inline V3      operator*(const double v) const; // Product of a number to a vector
			inline V3      operator/(const double v) const; // Division of a number to a vector
			inline V3      operator+(const V3 &v) const;    // Sum
			inline V3      operator-(const V3 &v) const;    // Substraction
			inline double  operator*(const V3 &v) const;    // Scalar product
			inline V3      operator^(const V3 &v) const;    // Vectorial product
			inline void    operator+=(const double v);      // Sum of a number to a vector
			inline void    operator-=(const double v);      // Substraction of a number to a vector
			inline void    operator*=(const double v);      // Product of a number to a vector
			inline void    operator/=(const double v);      // Division of a number to a vector
			inline void    operator+=(const V3 &v);         // Sum
			inline void    operator-=(const V3 &v);         // Substraction
			inline void    operator*=(const V3 &v);         // Point to point product

			inline bool    operator==(const V3 &v) const;   // Check equalty
			inline bool    operator<(const V3 &v) const;    // Check greater than

			inline double  operator[](int i) const;
			inline double &operator[](int i);

		private:
			double val[3];
	};

	// Constructor and destructor definitions
	inline V3::V3() { val[0] = val[1] = val[2] = 0.0; }
	inline V3::V3(const double x, const double y, const double z) { val[0]=x; val[1]=y; val[2]=z; }
	inline V3::V3(const double *x) { val[0]=x[0]; val[1]=x[1]; val[2]=x[2]; }
	inline V3::V3(const V3 &v) { val[0]=v.val[0]; val[1]=v.val[1]; val[2]=v.val[2]; }
	inline V3::~V3() { };

	// Operators
	inline V3     &V3::operator=(const V3 &v)             { val[0]=v[0]; val[1]=v[1]; val[2]=v[2]; return (*this); }
	inline V3      V3::operator+(const double v) const    { return( V3(val[0]+v,val[1]+v,val[2]+v) ); }
	inline V3          operator+(double val, const V3 &v) { return( V3(val+v[0],val+v[1],val+v[2]) ); }
	inline V3      V3::operator-(const double v) const    { return( V3(val[0]-v,val[1]-v,val[2]-v) ); }
	inline V3          operator-(double val, const V3 &v) { return( V3(val-v[0],val-v[1],val-v[2]) ); }
	inline V3      V3::operator*(const double v) const    { return( V3(val[0]*v,val[1]*v,val[2]*v) ); }
	inline V3          operator*(double val, const V3 &v) { return( V3(val*v[0],val*v[1],val*v[2]) ); }
	inline V3      V3::operator/(const double v) const    { return( V3(val[0]/v,val[1]/v,val[2]/v) ); }
	inline V3      V3::operator+(const V3 &v) const       { return( V3(val[0]+v[0],val[1]+v[1],val[2]+v[2]) ); }
	inline V3      V3::operator-(const V3 &v) const       { return( V3(val[0]-v[0],val[1]-v[1],val[2]-v[2]) ); }
	inline double  V3::operator*(const V3 &v) const       { return( val[0]*v[0] + val[1]*v[1] + val[2]*v[2] ); }
	inline V3      V3::operator^(const V3 &v) const       { return( V3(val[1]*v[2]-val[2]*v[1], -val[0]*v[2]+val[2]*v[0], val[0]*v[1]-val[1]*v[0])); }
	inline void    V3::operator+=(const double v)         { val[0]+=v; val[1]+=v; val[2]+=v; }
	inline void    V3::operator-=(const double v)         { val[0]-=v; val[1]-=v; val[2]-=v; }
	inline void    V3::operator*=(const double v)         { val[0]*=v; val[1]*=v; val[2]*=v; }
	inline void    V3::operator/=(const double v)         { val[0]/=v; val[1]/=v; val[2]/=v; }
	inline void    V3::operator+=(const V3 &v)            { val[0]+=v[0]; val[1]+=v[1]; val[2]+=v[2]; }
	inline void    V3::operator-=(const V3 &v)            { val[0]-=v[0]; val[1]-=v[1]; val[2]-=v[2]; }
	inline void    V3::operator*=(const V3 &v)            { val[0]*=v[0]; val[1]*=v[1]; val[2]*=v[2]; }
	inline bool    V3::operator==(const V3 &v) const      { V3 d = (*this) - v; return( (d.norm2() < 1e-10) ? true : false ); }
	inline bool    V3::operator<(const V3 &v)const        { return(this->norm2() < v.norm2()); }
	inline double  V3::operator[](int i) const            { return val[i]; }
	inline double &V3::operator[](int i)                  { return val[i]; }

	// Functions
	inline double  V3::V3::norm2() const { return(val[0]*val[0] + val[1]*val[1] + val[2]*val[2]); }
	inline void    V3::V3::print()       { std::printf("| %f %f %f |",val[0],val[1],val[2]); }

	// V3v:  Vector of V3
	class V3v {
		public:
			// Constructors and destructors
			inline V3v();
			inline V3v(const int nn);
			inline V3v(const int nn, const float *val);
			inline V3v(const int nn, const double *val);
			inline V3v(const int nn, const V3 *val);
			inline V3v(const V3v &v);
			inline ~V3v();

			// Functions
			inline void size(const int nn);
			inline void size(const int nn, const V3 *val);
			inline void size(const int nn, const float *val);
			inline void size(const int nn, const double *val);

			inline int     len();
			inline void   clear();
			inline V3     *data();
			inline float  *tofloat();
			inline double *todouble();
			inline bool    isempty();

			// Operators
			inline V3v &operator=(const V3v &v);
			
			inline V3   operator[](int i) const;
			inline V3  &operator[](int i);

			// Iterator
			class iterator {
				public:
					using value_type        = V3;
					using difference_type   = std::ptrdiff_t;
					using pointer           = V3*;
					using reference         = V3&;
					using _category         = std::random_access_iterator_tag;

					inline iterator() : v(nullptr), i(0) {}
					inline iterator(V3v* vv, int ii) : v(vv), i(ii) {}

					inline       V3&       operator*()            { return (*v)[i]; }
					inline const V3&       operator*() const      { return (*v)[i]; }
					inline       double&   operator[](int j)      { return (*v)[i][j]; }
					inline const double&   operator[](int j)const { return (*v)[i][j]; }

					inline       iterator& operator++()           { ++i; return *this; }
					inline       iterator& operator--()           { --i; return *this; }
					inline       iterator  operator++(int)        { iterator r(*this); ++i; return r; }
					inline       iterator  operator--(int)        { iterator r(*this); --i; return r; }

					inline       iterator& operator+=(int n)      { i += n; return *this; }
					inline       iterator& operator-=(int n)      { i -= n; return *this; }

					inline       iterator  operator+(int n) const { iterator r(*this); return r += n; }
					inline       iterator  operator-(int n) const { iterator r(*this); return r -= n; }

					inline difference_type operator-(iterator const& r) const { return i - r.i; }

					inline bool operator<(iterator const& r)  const { return i <  r.i; }
					inline bool operator<=(iterator const& r) const { return i <= r.i; }
					inline bool operator>(iterator const& r)  const { return i >  r.i; }
					inline bool operator>=(iterator const& r) const { return i >= r.i; }
					inline bool operator!=(const iterator &r) const { return i != r.i; }
					inline bool operator==(const iterator &r) const { return i == r.i; }

					inline int  ind() { return i; }

				private:
					V3v* v;
					int  i;		
			};

			inline iterator begin()                               { return V3v::iterator{this,0}; }
			inline iterator end()                                 { return V3v::iterator{this,n}; }

		private:
			int  n;
			V3  *v;
			bool alloc;
	};

	// Constructor and destructor definitions
	inline V3v::V3v()                                      { alloc = false; n = 0; }
	inline V3v::V3v(const int nn)                          { alloc = false; n = 0; size(nn); }
	inline V3v::V3v(const int nn, const float *val)        { alloc = false; n = 0; size(nn,val); }
	inline V3v::V3v(const int nn, const double *val)       { alloc = false; n = 0; size(nn,val); }
	inline V3v::V3v(const int nn, const V3 *val)           { alloc = false; n = 0; size(nn,val); }
	inline V3v::V3v(const V3v &val)                        { alloc = false; n = 0; size(val.n,val.v); }
	inline V3v::~V3v()                                     { clear(); }

	// Operators
	inline V3   V3v::operator[](int i) const               { return (i>=0) ? v[i] : v[n+i]; }
	inline V3  &V3v::operator[](int i)                     { return (i>=0) ? v[i] : v[n+i]; }

	inline V3v &V3v::operator=(const V3v &v)               { size(v.n,v.v); return (*this); }

	// Functions
	inline void V3v::size(const int nn)                    { n = nn; v = new V3[n]; alloc = true; }
	inline void V3v::size(const int nn, const V3 *val)     { size(nn); std::memcpy(v,val,n*sizeof(V3)); }
	inline void V3v::size(const int nn, const float *val)  { 
		size(nn);
		#pragma omp parallel for
		for(int i=0;i<n;i++)
			v[i] = V3((double)(val[3*i + 0]),(double)(val[3*i + 1]),(double)(val[3*i + 2]));
	}
	inline void V3v::size(const int nn, const double *val) { 
		size(nn);
		#pragma omp parallel for
		for(int i=0;i<n;i++)
			v[i] = V3(val[3*i + 0],val[3*i + 1],val[3*i + 2]);
	}
	inline int     V3v::len()                              { return n; }
	inline void    V3v::clear()                            { n = 0; if (alloc) { delete [] v; } alloc = false; }
	inline V3     *V3v::data()                             { return v; }
	inline bool    V3v::isempty()						   { return n == 0; }
	inline float  *V3v::tofloat()                          {
		float *out; out = new float[3*n];
		#pragma omp parallel for 
		for (int i=0;i<n;i++){
			out[3*i + 0] = (float)(v[i][0]);
			out[3*i + 1] = (float)(v[i][1]);
			out[3*i + 2] = (float)(v[i][2]);
		}
		return out;		
	}
	inline double *V3v::todouble()                         {
		double *out;  out = new double[3*n];
		#pragma omp parallel for
		for (int i=0;i<n;i++){
			out[3*i + 0] = v[i][0];
			out[3*i + 1] = v[i][1];
			out[3*i + 2] = v[i][2];
		}
		return out;		
	}
}

#endif
