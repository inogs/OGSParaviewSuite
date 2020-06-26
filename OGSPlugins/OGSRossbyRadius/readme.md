## OGS Rossby Radius

The **OGS Rossby Radius** computes the [Rossby radius of deformation](https://en.wikipedia.org/wiki/Rossby_radius_of_deformation), i.e., the length scale at which rotational effects become as important as buoyancy or gravity wave effects. This filter needs the seawater density field, which can be computed from the **OGS Water Density** (recommeded to work with potential density), and the mixing layer depth from **OGS Mixing Layer Depth**.

The barotropic (external) radius of deformation is defined as

<img src="https://github.com/inogs/OGSParaviewSuite/blob/master/OGSPlugins/OGSRossbyRadius/doc/RReq1.png" alt="" width="250"/>

where g=9.81 m/s2 is the gravitational acceleration, D is the maximum depth (bathymetry) and f=1e-4 1/s is the Coriolis frequency defined as

<img src="https://github.com/inogs/OGSParaviewSuite/blob/master/OGSPlugins/OGSRossbyRadius/doc/RReq2.png" alt="" width="250"/>

where Omega=7.2921e-5 rad/s is the rotation of the Earth and phi the latitude.

The internal Rossby radius of deformation is defined as

<img src="https://github.com/inogs/OGSParaviewSuite/blob/master/OGSPlugins/OGSRossbyRadius/doc/RReq3.png" alt="" width="250"/>

where MLD is the mixing layer depth ang g' the reduced gravity, defined as:

<img src="https://github.com/inogs/OGSParaviewSuite/blob/master/OGSPlugins/OGSRossbyRadius/doc/RReq4.png" alt="" width="250"/>

where rho1 is the averaged density in all the water column and rho2 is the averaged density on the MLD.