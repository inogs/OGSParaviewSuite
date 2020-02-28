/* 
	VORTEX
	
	Vortex class definition.

	Arnau Miro (OGS) (c) 2020
*/

#ifndef VORTEX_H
#define VORTEX_H

#include <cstdlib>
#include <cstring>
#include <vector>

#include "V3.h"
#include "field.h"

namespace vortex
{
	class Vortex {
		public:
			inline Vortex()                              { alloc = false; n = 0; }
			inline Vortex(const int nn, const int *e)    { alloc = false; set(nn); set_elems(e); }
			inline Vortex(const Vortex &v)               { alloc = false; set(v.n); set_elems(v.elms); }
			inline ~Vortex()                             { clear(); }

			inline void    set(const int nn)             { n = nn; allocElems(); alloc = true; }
			inline void    set_elems(const int *e)       { if (alloc) std::memcpy(elms,e,n*sizeof(int)); }
			inline void    clear()                       { deallocElems(); alloc = false; }
			inline int    *elems()                       { return elms; }
			inline int     get_n()                       { return n; }
			inline bool    isempty()                     { return n == 0; }

			inline v3::V3  baricenter(v3::V3v &xyz);
			inline double  size(v3::V3v &xyz);
			template <class T>
			inline int     icenter(v3::V3v &xyz, field::Field<T> &f);
			template <class T>
			inline v3::V3  center(v3::V3v &xyz, field::Field<T> &f);
			template <class T>
			inline v3::V3  rotation(field::Field<T> &f);
			template <class T>
			inline v3::V3  rotation(v3::V3v &xyz, field::Field<T> &f1, field::Field<T> &f2);
			template <class T>
			inline T       absolute_strength(field::Field<T> &w, field::Field<T> &f);
			template <class T>
			inline T       absolute_strength(v3::V3v &xyz, field::Field<T> &f1, field::Field<T> &f2);
			template <class T>
			inline T       relative_strength(field::Field<T> &w, field::Field<T> &f);
			template <class T>
			inline T       relative_strength(v3::V3v &xyz, field::Field<T> &f);

			inline int     operator[](int i) const       { return (i>=0) ? elms[i] : elms[n+i]; }
			inline int    &operator[](int i)             { return (i>=0) ? elms[i] : elms[n+i]; }
			inline Vortex &operator=(const Vortex &v)    { clear(); set(v.n); set_elems(v.elms); return (*this); }

		private:
			bool alloc;
			int n;
			int *elms;

			inline void allocElems()                     { if (!alloc) elms = new int[n]; }
			inline void deallocElems()                   { if (alloc)  delete [] elms; }
	};

	typedef std::vector<vortex::Vortex> VortexList;

	inline v3::V3 Vortex::baricenter(v3::V3v &xyz) {
		/*
			Compute the baricenter as the average of the
			points of the vortex.
		*/
		v3::V3 center(0.,0.,0.);
		// Loop the elements of the vortex
		for (int ii=0; ii<n; ++ii)
			center += xyz[elms[ii]];
		center /= (double)(n);
		return center;
	}

	inline double Vortex::size(v3::V3v &xyz) {
		v3::V3 center = this->baricenter(xyz);
		// Use a 2D approximation for the distance from
		// the center to the farthest extreme
		double dist = 0;
		for (int ii=0; ii<n; ++ii) {
			v3::V3 d  = xyz[elms[ii]] - center; d[2] = 0;
			double dd = d.norm2();
			dist = (dd > dist) ? dd : dist;
		}
		return sqrt(dist);
	}

	template <class T>
	inline int Vortex::icenter(v3::V3v &xyz, field::Field<T> &f) {
		// Compute the center as the position where the maximum
		// value of the field "f" is found within the vortex.
		// Further refinements might apply.
		T valmax = -1e-30; int posmax = -1;
		for (int ii=0; ii<n; ++ii) {
			if (f[elms[ii]][0] > valmax) { 
				valmax = f[elms[ii]][0]; 
				posmax = elms[ii];
			}
		}
		return posmax;
	}

	template <class T>
	inline v3::V3 Vortex::center(v3::V3v &xyz, field::Field<T> &f) {
		// Compute the center as the position where the maximum
		// value of the field "f" is found within the vortex.
		// Further refinements might apply.
		return xyz[this->icenter<T>(xyz,f)];
	}

	template <class T>
	inline v3::V3 Vortex::rotation(field::Field<T> &f) {
		// Compute the orientation of the rotation. Returns a normalized
		// vector indicative of the rotation.
		// f must be a vector
		v3::V3 rot(0.,0.,0.);
		// Loop the elements of the vortex
		for (int ii=0; ii<n; ++ii)
			rot += v3::V3(f[elms[ii]][0],f[elms[ii]][1],f[elms[ii]][2]);
		rot /= sqrt(rot.norm2());
		return rot;
	}

	template <class T>
	inline v3::V3 Vortex::rotation(v3::V3v &xyz, field::Field<T> &f1, field::Field<T> &f2) {
		// Compute the orientation of the rotation. Returns a normalized
		// vector indicative of the rotation. Takes the value of the center
		// f must be a vector
		int ii = this->icenter<T>(xyz,f1);
		v3::V3 rot(f2[ii][0],f2[ii][1],f2[ii][2]);
		rot /= sqrt(rot.norm2());
		return rot;
	}

	template <class T>
	inline T Vortex::absolute_strength(field::Field<T> &w, field::Field<T> &f) {
		// Compute the absolute strength of the vortex by
		// performing a weighted average on the field f
		T sum_weight = 0., meanval = 0.;
		for (int ii=0; ii<n; ++ii) {
			T v = sqrt(f[elms[ii]][0]*f[elms[ii]][0] + f[elms[ii]][1]*f[elms[ii]][1] 
				+ f[elms[ii]][2]*f[elms[ii]][2]);
			sum_weight += w[elms[ii]][0];
			meanval    += w[elms[ii]][0]/sum_weight*(v-meanval);
		}
		return meanval;
	}

	template <class T>
	inline T Vortex::absolute_strength(v3::V3v &xyz, field::Field<T> &f1, field::Field<T> &f2) {
		// Compute the absolute strength of the vortex by
		// returning the value at the center
		int ii = this->icenter<T>(xyz,f1);
		T val  = sqrt(f2[ii][0]*f2[ii][0] + f2[ii][1]*f2[ii][1] + f2[ii][2]*f2[ii][2]);
		return val;
	}

	template <class T>
	inline T Vortex::relative_strength(field::Field<T> &w, field::Field<T> &f) {
		// Compute the relative strength of the vortex by
		// performing a weighted average on the field f
		T sum_weight = 0., meanval = 0.;

		for (int ii=0; ii<n; ++ii) {
			sum_weight += w[elms[ii]][0];
			meanval    += w[elms[ii]][0]/sum_weight*(f[elms[ii]][0]-meanval);
		}
		return meanval;
	}	

	template <class T>
	inline T Vortex::relative_strength(v3::V3v &xyz, field::Field<T> &f) {
		// Compute the relative strength of the vortex by
		// returning the value at the center
		int ii = this->icenter<T>(xyz,f);
		return f[ii][0];
	}
}

#endif