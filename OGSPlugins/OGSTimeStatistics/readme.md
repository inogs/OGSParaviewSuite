## OGS Time Statistics

The **OGS Time Statistics** are a set of plugins that perform temporal averages of a dataset. They share an averaging algorithm based on the parallel 1 step algorithm for the computation of the mean (see **OGS Spatial Stats**). There are two modes of of operation:
1. A default versatile mode that interates the ParaView pipeline through the requested instants to perform the average, thus at each instant the pipeline needs to be recomputed. This algorithm is parallelized using OpenMP.
2. A parallel MPI algorithm that loads directly the NetCDF files, hence, it must be applied after the **OGS Reader** and needs a _vtkRectilinearGrid_.

### OGS Time Average

The **OGS Time Average** performs the temporal average of a dataset given an interval of time. With that, the plugin uses a generic time requestor to search the data time list (generated with the time-stamps) and obtain which data must be averaged and their weights.

The plugin writes just one time-step, therefore, all the temporal variability is lost. Nevertheless, the plugin does not destroy the time information so it can be used in conjunction with other time-dependent plugins.

### OGS Time Aggregator

The **OGS Time Aggregator** allows to convert data from one time frequency to another (e.g., monthly to yearly). This plugin uses the TimeList method getXXXList (where XXX stands for yearly, monthly, etc.) to generate a list of requestors from where to extract the instants and weights for the average.

The plugin generates a new timestep list based on the requestors list. The new ParaView timestep controls which average needs to be performed, then, according to this instant, the list of instants and weights is obtained. This means that the averages are performed "on demand" every time the user moves through time.

### OGS Climatology

The **OGS Climatology** allows to aggregate the data climatologically, i.e., group it by season, month or day regardless of the year. This plugin uses the climatological requestors to generate a list of instants and weights from where to perform the averages.

The plugin generates new timesteps based on the climatological analysis requested (4 for seasons, 12 for months and 365 for days). Similar to the **OGS Time Aggregator**, the new timesteps control which instants are selected for averaging, hence the averages are done "on demand" every time the user moves through time.