## OGS Reader

The **OGS Reader** is a ParaView plugin that lets import data from [OGSTM-BFM](http://medeaf.inogs.it/how-data-are-generated) and/or [Copernicus Marine](http://marine.copernicus.eu/services-portfolio/access-to-products/). The plugin needs to read information from a masterfile (_.ogs_), a mesh file (_.ogsmsh_) and NetCDF files (_.nc_) containing the variables.

### Generation of the masterfile

Generating the masterfile and the meshes is done by means of the **OGS2Paraview** script, which is a helper python script included in the suite and deployed during the installation. It can be called by issuing:

```bash
OGS2ParaView [-h] -n NAME -p PATH [-m MODE] [-c CONF] [--gen-mesh] [--meshmask MESH]
```
where _NAME_ is the name of the masterfile and _PATH_ is the path to the dataset. Four _MODE_ are available:

* 0: Generate the mesh and exit
* 1: Use AVE_PHYS to generate the timesteps
* 2: Use AVE_FREQ_1 to generate the timesteps
* 3: Use FORCINGS to generate the timesteps
* 4: Use GENERALS to generate the timesteps

A recommended directory tree for the simulation would be:

* _meshmask.nc_
* _AVE_PHYS_: directory containing the physical variables (ave.\*.phys.nc)
* _AVE_FREQ_: directory containing the biogeochemical variables (ave.\*.\*.nc)
* _FORCINGS_: directory containing the simulation forcings
* _GENERALS_: directory containing any kind of variable.

The following rules are recommended to be observed when organizing the data directories:

* The variables must be of the same dimensions specified in the meshmask.
* There should be all the variables for a specific timestep.
* Variable file names should be consistent (<prefix>.<date>.<name>.nc).

The data directory and the variable file names can be modified by specifying a configuration file when running the tool.

### Generating the mesh files

The first consideration when generating the mesh is to project the angular coordinates (longitude and latitude) to meters. To reconstruct the mesh, two coordinates are needed:

* Coordinates of the cell vertices, used to generate the _vtkRectilinearGrid_ (glamf, gphif, gdepw).
* Coordinates of the cell centers, where to represent the variables (glamt,gphit,nav_lev)

This projection is done using the basemap library from Python. By default the Mercator 

<img src="https://github.com/inogs/OGSParaviewSuite/blob/master/OGSPlugins/OGSReader/doc/eq1.png" alt="" width="500"/>

with

<img src="https://github.com/inogs/OGSParaviewSuite/blob/master/OGSPlugins/OGSReader/doc/eq2.png" alt="" width="200"/>

and the cylindrical 

<img src="https://github.com/inogs/OGSParaviewSuite/blob/master/OGSPlugins/OGSReader/doc/eq3.png" alt="" width="180"/>

projections are considered.

<img src="https://github.com/inogs/OGSParaviewSuite/blob/master/OGSPlugins/OGSReader/doc/mesh.png" alt="" width="600"/>

The basins mask and the continental shelf mask (coast mask) are also included in the mesh file. The basins mask is a mask that subdivides the different regions of the Mediterranean and follows the battimetry. The continental shelf mask subdivides the regions of the Mediterranean whose depth is less than 200 meters.

In addition to this mask arrays, tensor fields for e1, e2 and e3 are also created. These contain information on the mesh stretching and are necessary to project magnitudes to the cell centers, compute gradients or compute spatial statistics. They are stored as:

| eij |  0  |  1  |  2  |  3  |
| --- | --- | --- | --- | --- |
|  1  | e1t | e1u | e1v | e1f |
|  2  | e2t | e2u | e2v | e2f |
|  3  | e3t | e3u | e3v | e3w |

### Importing into ParaView

The variables inside the NetCDF files are read and stored as _fields_. Then, they are converted into _vtk_ using the _vtkField_ interface. Special care must be taken with the velocity as it is defined at the cell faces and must be projected to the cell centers. Moreover, the depth levels are multiplied by a scaling factor, otherwise, the data would seem two dimensional.

The _Metadata_ string array contains information useful to the visualization. It is important to carry the _Metadata_ in each operation, as it is usually used to recover the scaling factor or the timestep.