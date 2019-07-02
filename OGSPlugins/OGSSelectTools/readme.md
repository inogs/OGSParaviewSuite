## OGS Select Tools

The **OGS Select Tools** plugin lets the user select between several sub-basins or between the continental shelf and the open sea. 

The code generates a cut mask and populates it as _true_ for the values of the mesh that the user has selected. Then, this cut mask is feeded to the _vtkThreshold_ algorithm; which cuts all values from 0.5 to 1 (i.e., higher than 0). This filter can operate with either cell or point data.

A particularity of this filter is that it returns a _vtkUnstructuredGrid_, thus further forbidding any operations on a _vtkRectilinearGrid_.

The masks use the _uint8_t_ type for an efficient storage of small integer arrays. In particular, in order to avoid overlapping issues, the basins mask is formed by a tensor of 16 components being each one the representation of a basin. The magnitude of the tensor potentially shows the Mediterranean sea with the reginos of overlapping highlighted.
