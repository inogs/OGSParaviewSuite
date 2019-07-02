## OGS Plot View

The **OGS Plot Views** are a collection of python plot views that perform specific plots that are not currently programmed inside ParaView. They are based on a _vtkPythonView_ with specific Python scripts and controls based on the _vtkPythonProgrammableFilter_. These views are:

* **OGS Vertical Profile Plot**: this view is used with the **OGS Depth Profile** to provide a plot of the depth profile with the depth on the y axis and the correct vertical scaling.
* **OGS Hovmoeller Plot**: this view is used with the **OGS Hovmoeller** and enables to produce Hovmoeller plots using matplotlib.
* **OGS Spaghetti Plot**: this view is used with the **OGS Spaghetti** to provide a spaghetti plot with the time at the _x axis_,

It is also possible to produce a high qualitiy PFD from this views.