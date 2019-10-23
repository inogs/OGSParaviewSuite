/* 
    FIELD
	
	Definition of a generic field for CFD usage.
	Contains some basic operations and easy access.

	Arnau Miro (UPC-ESEIAAT) (c) 2018
*/

#ifndef Field_h
#define Field_h

#include <cstdlib>
#include <cstring>
#include <algorithm>

namespace field
{
	template <class T>
	class Field;

	template<class T> 
	Field<T> operator+(const T v, const Field<T> &ff);
	template<class T> 
	Field<T> operator-(const T v, const Field<T> &ff);
	template<class T> 
	Field<T> operator*(const T v, const Field<T> &ff);
	template<class T> 
	Field<T> operator/(const T v, const Field<T> &ff);

	template <class T>
	class Field {

		friend Field<T> operator+<>(const T v, const Field<T> &ff);
		friend Field<T> operator-<>(const T v, const Field<T> &ff);
		friend Field<T> operator*<>(const T v, const Field<T> &ff);
		friend Field<T> operator/<>(const T v, const Field<T> &ff);

		public:
			// Constructors and destructors
			inline Field()                                            { alloc = false; n = 0; m = 0; sz = 0; }
			inline Field(const int nn, const int mm)                  { alloc = false; set_dim(nn,mm); }
			inline Field(const int nn, const int mm, const T  v)      { alloc = false; set(nn,mm,v); }
			inline Field(const int nn, const int mm, const T *v)      { alloc = false; set(nn,mm,v); }
			inline Field(const Field<T> &f)                           { alloc = false; set(f.n,f.m,f.val); }
			inline ~Field()                                           { clear(); }

			// Functions
			inline void   set_dim(const int nn, const int mm)         { n = nn; m = mm; sz = (size_t)(n*m); val = new T[sz]; alloc = true; }
			inline void   set_val(const T  v)                         { if (alloc) std::fill(val,val+sz,v); }
			inline void   set_val(const T *v)                         { if (alloc) std::memcpy(val,v,sz*sizeof(T)); }
			inline void   clear()                                     { n = 0; if (alloc) { delete [] val; } alloc = false; }
			inline void   set(const int nn, const int mm, const T  v) { set_dim(nn,mm); set_val(v); }
			inline void   set(const int nn, const int mm, const T *v) { set_dim(nn,mm); set_val(v); }
			inline T     *data()                                      { return val; }
			inline int    get_n()                                     { return n; }
			inline int    get_m()                                     { return m; }
			inline size_t get_sz()                                    { return sz; }
			inline bool   isempty()                                   { return n == 0; }

			template <class U>
			inline Field<U> convert()                                 { Field<U> f(n,m); std::copy(val,val+sz,f.data()); return f; }

			// Operators
			inline T         *operator[](int i)                   { return (i>=0) ? val + m*i : val + m*(n+i); }
			inline Field<T>  &operator=(const Field<T> &f)        { set(f.n,f.m,f.val); return (*this); }

			inline Field<T>   operator+(const T v) const          { Field<T> f(n,m); for(int i=0;i<sz;i++) {f.val[i] = val[i] + v;} return f; }
			inline Field<T>   operator-(const T v) const          { Field<T> f(n,m); for(int i=0;i<sz;i++) {f.val[i] = val[i] - v;} return f; }
			inline Field<T>   operator*(const T v) const          { Field<T> f(n,m); for(int i=0;i<sz;i++) {f.val[i] = val[i] * v;} return f; }
			inline Field<T>   operator/(const T v) const          { Field<T> f(n,m); for(int i=0;i<sz;i++) {f.val[i] = val[i] / v;} return f; }
			inline Field<T>   operator+(const Field<T> &ff) const { Field<T> f(n,m); for(int i=0;i<sz;i++) {f.val[i] = val[i] + ff.val[i];} return f; }
			inline Field<T>   operator-(const Field<T> &ff) const { Field<T> f(n,m); for(int i=0;i<sz;i++) {f.val[i] = val[i] - ff.val[i];} return f; }
			inline Field<T>   operator*(const Field<T> &ff) const { Field<T> f(n,m); for(int i=0;i<sz;i++) {f.val[i] = val[i] * ff.val[i];} return f; }
			inline Field<T>   operator/(const Field<T> &ff) const { Field<T> f(n,m); for(int i=0;i<sz;i++) {f.val[i] = val[i] / ff.val[i];} return f; }
			inline void       operator+=(const T v)               { for(int i=0;i<sz;i++) {val[i] += v;} }
			inline void       operator-=(const T v)               { for(int i=0;i<sz;i++) {val[i] -= v;} }
			inline void       operator*=(const T v)               { for(int i=0;i<sz;i++) {val[i] *= v;} }
			inline void       operator/=(const T v)               { for(int i=0;i<sz;i++) {val[i] /= v;} }
			inline void       operator+=(const Field<T> &f)       { for(int i=0;i<sz;i++) {val[i] += f.val[i];} }
			inline void       operator-=(const Field<T> &f)       { for(int i=0;i<sz;i++) {val[i] -= f.val[i];} }
			inline void       operator*=(const Field<T> &f)       { for(int i=0;i<sz;i++) {val[i] *= f.val[i];} }
			inline void       operator/=(const Field<T> &f)       { for(int i=0;i<sz;i++) {val[i] /= f.val[i];} }

			inline bool       operator==(const Field<T> &f) const { return(n == f.n); }
			inline bool       operator!=(const Field<T> &f) const { return(n != f.n); }
			inline bool       operator<(const Field<T> &f) const  { return(n < f.n); }
			inline bool       operator>(const Field<T> &f) const  { return(n > f.n); }

			// Iterator
			class iterator {
				public:
					using value_type        = T;
					using difference_type   = std::ptrdiff_t;
					using pointer           = T*;
					using reference         = T&;
					using _category = std::random_access_iterator_tag;

					inline iterator() : f(nullptr), i(0) {}
					inline iterator(Field<T>* ff, int ii) : f(ff), i(ii) {}

					inline       pointer   operator*()            { return (*f)[i]; }
					inline const pointer   operator*() const      { return (*f)[i]; }
					inline       reference operator[](int j)      { return (*f)[i][j]; }
					inline       reference operator[](int j)const { return (*f)[i][j]; }

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
					Field<T>* f;
					int       i;		
			};

			inline iterator begin()                               { return Field<T>::iterator{this,0}; }
			inline iterator end()                                 { return Field<T>::iterator{this,n}; }

		private:
			int     n;     // length of the field
			int     m;     // number of elements of the field
			T      *val;   // values as a 1D C array
			size_t  sz;    // length of the vector (sz = n*m)
			bool    alloc; // has it been allocated?
	};

	template<class T> 
	Field<T> operator+(const T v, const Field<T> &ff) { Field<T> f(ff.n,ff.m); for(int i=0;i<ff.sz;i++) {f.val[i] = v + ff.val[i];} return f; }
	template<class T> 
	Field<T> operator-(const T v, const Field<T> &ff) { Field<T> f(ff.n,ff.m); for(int i=0;i<ff.sz;i++) {f.val[i] = v - ff.val[i];} return f; }
	template<class T> 
	Field<T> operator*(const T v, const Field<T> &ff) { Field<T> f(ff.n,ff.m); for(int i=0;i<ff.sz;i++) {f.val[i] = v / ff.val[i];} return f; }
	template<class T> 
	Field<T> operator/(const T v, const Field<T> &ff) { Field<T> f(ff.n,ff.m); for(int i=0;i<ff.sz;i++) {f.val[i] = v * ff.val[i];} return f; }
}

#endif 