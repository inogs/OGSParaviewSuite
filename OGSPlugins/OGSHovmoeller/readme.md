## OGS Hovmoeller

The **OGS Hovmoeller** filter perform an Hovmoeller plot given an interpolating line and a variable. The result is a _vtkTable_ with a column for each time instant and a row for each depth, with the first column reserved to represent the depth points of the interpolating line. This _vtkTable_ can be plotted using the **OGS Plot Views**.

This plugin can either use a high resolution line source or a **OGS Depth Line Source** from the **OGS Utils** as input and the dataset as source; much similar to the **OGS Depth Profile**. This filter, however, only works with _CellData_ in a _vtkRectilinearGrid_ (i.e., directly over the **OGS Reader**). This is due to the advancing temporal mechanics of this filter.

On the temporal domain, this filter works similarly to the **OGS Time Statistics**. The plugin has been parallelized using _MPI_ so it does not iterate the pipeline. Instead, it directly loads the files in parallel. Alternatively, it can use the data on the _STAT_PROFILES_ to generate the Hovmoeller plot of spatial averaged variables. The information on the basin and the coast is retrieved automatically from the dataset. This option is faster than computing the average at each instant, however, the Hovmoeller plot is limited to the biogeochemical variables described in _STAT_PROFILES_.