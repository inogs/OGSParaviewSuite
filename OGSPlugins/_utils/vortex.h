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
			inline Vortex()                              { alloc = false; ic = -1; n = 0; }
			inline Vortex(const int nn, const int *e)    { alloc = false; ic = -1; set(nn); set_elems(e); }
			inline Vortex(const Vortex &v)               { alloc = false; ic = -1; set(v.n); set_elems(v.elms); }
			inline ~Vortex()                             { clear(); }

			inline void    set(const int nn)             { n = nn; allocElems(); alloc = true; }
			inline void    set_elems(const int *e)       { if (alloc) std::memcpy(elms,e,n*sizeof(int)); }
			inline void    clear()                       { deallocElems(); alloc = false; }
			inline int    *elems()                       { return elms; }
			inline int     get_n()                       { return n; }
			inline bool    isempty()                     { return n == 0; }

			inline v3::V3 baricenter(v3::V3v &xyz);
			template <class T>
			inline int    icenter(field::Field<T> &f);
			template <class T>
			inline v3::V3 center(v3::V3v &xyz, field::Field<T> &f);
			template <class T>
			inline T      size(v3::V3v &xyz);
			template <class T>
			inline T      size(v3::V3v &xyz, field::Field<T> &f);
			template <class T>
			inline v3::V3 rotation(field::Field<T> &f);
			template <class T>
			inline T      absolute_strength(field::Field<T> &f);
			template <class T>
			inline T      absolute_strength(field::Field<T> &w, field::Field<T> &f);
			template <class T>
			inline T      relative_strength(field::Field<T> &f);
			template <class T>
			inline T      relative_strength(field::Field<T> &w, field::Field<T> &f);
			template <class T>
			inline T      circulation(field::Field<T> &no, field::Field<T> &f);

			inline int     operator[](int i) const       { return (i>=0) ? elms[i] : elms[n+i]; }
			inline int    &operator[](int i)             { return (i>=0) ? elms[i] : elms[n+i]; }
			inline Vortex &operator=(const Vortex &v)    { clear(); set(v.n); set_elems(v.elms); return (*this); }

		private:
			bool alloc;
			int n, ic;
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

	template <class T>
	inline int Vortex::icenter(field::Field<T> &f) {
		/*
			Compute the center as the position where the maximum
			value of the field "f" is found within the vortex.
			Further refinements might apply.
		*/
		T valmax = -1e-30; ic = -1;
		for (int ii=0; ii<n; ++ii) {
			if (f[elms[ii]][0] > valmax) { 
				valmax = f[elms[ii]][0]; 
				ic     = elms[ii];
			}
		}
		return ic;
	}

	template <class T>
	inline v3::V3 Vortex::center(v3::V3v &xyz, field::Field<T> &f) {
		/*
			Compute the center as the position where the maximum
			value of the field "f" is found within the vortex.
			Further refinements might apply.
		*/
		if (ic < 0) ic = this->icenter<T>(f);
		return xyz[ic];
	}

	template <class T>
	inline T Vortex::size(v3::V3v &xyz) {
		/*
			Compute the size as the furthest distance to the center.
		*/
		v3::V3 center = (ic<0) ? this->baricenter(xyz) : xyz[ic];
		T dist = 0.;
		for (int ii=0; ii<n; ++ii) {
			v3::V3 d  = (xyz[elms[ii]] - center);
			T dd = (T)(d.norm2());
			dist = (dd > dist) ? dd : dist;
		}
		return sqrt(dist);
	}

	template <class T>
	inline T Vortex::size(v3::V3v &xyz, field::Field<T> &f) {
		/*
			Compute the size as the distance to where the
			relative strength has decreased by 95%.
		*/
		v3::V3 center = (ic<0) ? this->center(xyz,f) : xyz[ic];
		T flim = (1-0.95)*f[ic][0], dist = -1e30;
		for (int ii=0; ii<n; ++ii) {
			if (f[elms[ii]][0] < flim) continue;
			v3::V3 d  = (xyz[elms[ii]] - center);
			T dd = (T)(d.norm2());
			dist = (dd > dist) ? dd : dist;
		}
		return sqrt(dist);
	}

	template <class T>
	inline v3::V3 Vortex::rotation(field::Field<T> &f) {
		/*
			Compute the orientation of the rotation. Returns a normalized
			vector indicative of the rotation. f must be a vector array.
		*/
		v3::V3 rot(0.,0.,0.);
		if (ic<0) {
			// Loop the elements of the vortex
			for (int ii=0; ii<n; ++ii)
				rot += v3::V3(f[elms[ii]][0],f[elms[ii]][1],f[elms[ii]][2]);
		} else {
			rot = v3::V3(f[ic][0],f[ic][1],f[ic][2]);
		}
		rot /= sqrt(rot.norm2());
		return rot;
	}

	template <class T>
	inline T Vortex::absolute_strength(field::Field<T> &f) {
		/*
			Compute the absolute strength of the vortex by 
			returning the value at the center.
		*/
		if (ic < 0) return 0.;
		T val = (f.get_m() == 1) ? f[ic][0] : sqrt(f[ic][0]*f[ic][0] + f[ic][1]*f[ic][1] + f[ic][2]*f[ic][2]);
		return val;
	}

	template <class T>
	inline T Vortex::absolute_strength(field::Field<T> &w, field::Field<T> &f) {
		/*
			Compute the absolute strength of the vortex by
			performing a weighted average on the field f.
		*/
		T sum_weight = 0., meanval = 0.;
		for (int ii=0; ii<n; ++ii) {
			T v = (f.get_m() == 1) ? f[elms[ii]][0] : 
				sqrt(f[elms[ii]][0]*f[elms[ii]][0] + f[elms[ii]][1]*f[elms[ii]][1] + f[elms[ii]][2]*f[elms[ii]][2]);
			sum_weight += w[elms[ii]][0];
			meanval    += w[elms[ii]][0]/sum_weight*(v-meanval);
		}
		return meanval;
	}

	template <class T>
	inline T Vortex::relative_strength(field::Field<T> &f) {
		/*
			Compute the relative strength of the vortex by
			returning the value at the center.
		*/
		if (ic < 0) return 0.;
		return f[ic][0];
	}

	template <class T>
	inline T Vortex::relative_strength(field::Field<T> &w, field::Field<T> &f) {
		/*
			Compute the relative strength of the vortex by
			performing a weighted average on the field f.
		*/
		T sum_weight = 0., meanval = 0.;

		for (int ii=0; ii<n; ++ii) {
			sum_weight += w[elms[ii]][0];
			meanval    += w[elms[ii]][0]/sum_weight*(f[elms[ii]][0]-meanval);
		}
		return meanval;
	}

	template <class T>
	inline T Vortex::circulation(field::Field<T> &no, field::Field<T> &f) {
		/*
			Compute the circulation of the vortex by 
			performing the integral of -f*dA.
		*/
		T Gamma = 0.;
		for (int ii=0; ii<n; ++ii) {
			int iel = elms[ii];
			// Scalar product of f*dA
			Gamma += (no[iel][0]*f[iel][0] + no[iel][1]*f[iel][1] + no[iel][2]*f[iel][2]);
		}
		return Gamma;
	}
}

#endif
