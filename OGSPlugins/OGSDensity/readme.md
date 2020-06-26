## OGS Water Density

The **OGS Water Density** implements the computation of the potential sea water density (in kg/m3) from the temperature and salinity. Optionally, the user can activate the _in-situ_ option that also uses the depth (or pressure) in the computation.

The following methods are implemented:
* Linear relationship with the temperature (_does not allow in-situ_).
* Linear relationship with the temperature and salinity (_does not allow in-situ_).
* Jackett and McDougall 1994.
* JAOT 12, Jackett and McDougall 1995 modified UNESCO polynomial.
* JAOT 20, Jackett et al. 2003

To be implemented is the [TEOS-10](http://www.teos-10.org/), the Thermodynamic Equation Of Seawatert 2010 that has been adopted by the Intergovernmental Oceanographic Comission.