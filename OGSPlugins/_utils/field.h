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
			inline Field()                                          { alloc = false; }
			inline Field(const int nn, const int mm)                { alloc = false; set_dim(nn,mm); }
			inline Field(const int nn, const int mm, const T *v)    { alloc = false; set(nn,mm,v); }
			inline Field(const Field<T> &f)                         { alloc = false; set(f.n,f.m,f.val); }
			inline ~Field()                                         { if (alloc) delete [] val; }

			// Functions
			inline void   set_dim(const int nn, const int mm)         { n = nn; m = mm; sz = (size_t)(n*m); val = new T[sz]; alloc = true;}
			inline void   set_val(const T *v)                         { std::memcpy(val,v,sz*sizeof(T)); }
			inline void   set(const int nn, const int mm, const T *v) { set_dim(nn,mm); set_val(v); }
			inline T     *data()                                      { return val; }
			inline int    get_n()                                     { return n; }
			inline int    get_m()                                     { return m; }
			inline size_t get_sz()                                    { return sz; }
			inline bool   isempty()                                   { return val == NULL;}

			// Operators
			inline T         *operator[](int i)                   { return val + m*i; }
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