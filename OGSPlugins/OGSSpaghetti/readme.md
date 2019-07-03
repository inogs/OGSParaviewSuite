## OGS Spaghetti

The **OGS Spaghetti** plugin performs a spaghetti plot given a point on the dataset and a variable. The result is a _vtkTable_ with a column for each time instant.

The plugin uses a point source as input and works similarly to the **OGS Hovmoeller**. This filter only works for cell data in a _vktRectilinearGrid_, as it needs to load the variables from NetCDF files at each instant. On the temporal domain, this filter works in a similar manner to the **OGS Temporal Statistics**. The plugin has been parallelized using MPI and it does not iterate over the pipeline.

Alternatively, as the **OGS Hovmoeller**, this filter can use the data on the _STAT_PROFILES_ to generate plots with the spatial averaged variables. The information on the basin and the coast is retrieved automatically from the dataset. This option is faster than iterating the pipeline since the variables have already been packed temporally; however, it is limited to the biogeochemical variables.

This filter can also produce and average on depth by performing the mean of the point of an interpolating line. Since the points on the interpolating line are equally spaced, a weighted average is obtained when elements of different height are crossed. The same idea applies when the line finishes at a depth that is not the boundary between two elements.