## OGS Mixing Layer Depth

The **OGS Mixing Layer Depth** computes the mixed layer where active turbulence has homogenized some range of depths.

Two criterias have been implemented, one using the temperature and the other using density:
* Temperature: |T-Tref| > dTref
* Density: |rho-rho_ref| > drho_ref
where Tref and rho_ref are taken at z=10 m by default and dTref=0.2; drho_ref=0.03 (de Boyer Montégut, 2004).

For performance reasons, this filter should be applied to a _vtkRectilinearGrid_ before converting to a _vtkUnstructuredGrid_, i.e., before applying any of the ** OGS Select Tools**.