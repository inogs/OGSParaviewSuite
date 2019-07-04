## OGS Time Statistics

The **OGS Time Statistics** plugin performs the temporal average of a dataset given the iniital and final time frames. This plugin uses the time stamp in the time steps to generate the strings for the time interval selector. Moreover, it is also uses the parallel 1 step algorithm for the computation of the mean (see **OGS Spatial Stats**).

This plugin is parallelized using MPI; each rank deals with a part of the dataset. This plugin does not iterate over the pipeline, instead it loads directly the NetCDF files. Hence, this plugin must be applied on the **OGS Reader** and needs a _vtkRectilinearGrid_.

This algorithm consists of four phases:

1. Initialize; where the arrays are build and broadcasted if needed.
2. Accumulate; where the average is computed.
3. Reduction; only done if run in parallel.
4. Finalize; where the output is stored.