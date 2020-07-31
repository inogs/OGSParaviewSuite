## OGS Plot Views

The **OGS Plot Views** are a collection of python plot views that perform specific plots that are not currently programmed inside ParaView. They are based on a _vtkPythonView_ with specific Python scripts and controls based on the _vtkPythonProgrammableFilter_. These views are:

* **OGS Vertical Profile Plot**: this view is used with the **OGS Depth Profile** to provide a plot of the depth profile with the depth on the y axis and the correct vertical scaling.
* **OGS Hovmoeller Plot**: this view is used with the **OGS Hovmoeller** and enables to produce Hovmoeller plots using matplotlib.
* **OGS Spaghetti Plot**: this view is used with the **OGS Spaghetti** to provide a spaghetti plot with the time at the _x axis_ and the evolution of a variable on the _y axis_
* **OGS Map Plot**: this view can be used to visualize slice point data (using _Cell Data to Point Data_) on a map. It uses [**Cartopy**](https://scitools.org.uk/cartopy/docs/latest/) to handle the projection, so it must be installed on the system in order to use this plot.

It is also possible to produce a high qualitiy PNG and PFD from this views. The NASA Blue Marble map is included for the **OGS Map Plot** as a high resolution Earth map.