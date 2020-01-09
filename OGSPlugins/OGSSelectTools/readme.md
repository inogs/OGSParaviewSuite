## OGS Select Tools

The **OGS Select Tools** plugin lets the user select between several sub-basins, between the continental shelf and the open sea or the Okubo-Weiss regions. It also enables the user to select the land areas and perform a cut using a user selected polygon.

The code generates a cut mask and populates it as _true_ for the values of the mesh that the user has selected. Then, this cut mask is feeded to the _vtkThreshold_ algorithm; which cuts all values from 0.5 to 1 (i.e., higher than 0). This filter can operate with either cell or point data. The polygon for the **OGS Select Polygon** filter can be generated using the _poly_selector_ tool from _bit.sea_. A particularity of this filter is that it returns a _vtkUnstructuredGrid_, thus further forbidding any operations on a _vtkRectilinearGrid_.

The masks use the _uint8_t_ type for an efficient storage of small integer arrays. In particular, in order to avoid overlapping issues, the basins mask is formed by a tensor of 16 components being each one the representation of a basin. The magnitude of the tensor potentially shows the Mediterranean sea with the reginos of overlapping highlighted.

The **OGS Land Outline** is a composed filter defined as a **OGS Select Land**, a **Slice** to cut the surface and a **Feature Edges** filter to compute the silhouette. Since it is a composed filter, its installation is different and optional. In order to install it, go to _Tools/Manage Custom Filters..._ and _Import_ the XML file (_OGSLandOutline.xml_).