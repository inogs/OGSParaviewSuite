## OGS Spatial Stats

The **OGS Spatial Stats** filter is able to generate the following nine statistics to any part (or the whole) Mediterranean sea:

* Mean
* Standard deviation (std)
* Minimum (min)
* Maximum (max)
* Percentile 5% (p05)
* Percentile 25% (p25)
* Percentile 50% (p50)
* Percentile 75% (p75)
* Percentile 95% (p95)

Although this fitler can be used under any kind of data, cell data is strongly advised. 

By default, the weighting used is the cell area (e1te2t) although the user can select to use the cell volume (e1te2te3t). In any case, it is necessary to load the mesh weights for this filter to work.

Moreover, as default, the code generated an average per depth level (that is scanned at the beginning). Alternatively, the user can decide to provide depth levels by inputting them in the selection box. Note that the levels go from 0 to the selected value and so on until the bottom is reached.

This filter has  been validated against the **OGS Spatial Stats from File** filter and produces the same results under the same circumstances.

### Computation of the statistics

The computation of the statistics deserves a special mention. The general naive approach to compute the standard deviation takes two loops on the mesh: a first one for the mean and a second for the variance.

Welford (1962) proposed an algorithm to comptue the variance in a single loop. West (1979) extended the method for weighted computations using a recursive algorithm:

<img src="https://github.com/inogs/OGSParaviewSuite/blob/master/OGSPlugins/OGSSpatialStats/doc/eq1.png" alt="" width="200"/>

where the mean and standard deviation can be computed acummulating over _i_ until _n_.

A parallel algorithm was proposed by Chan (1979). the magnitudes are basically computed as described before for each sub-domain and the following recursive reduction is used:

<img src="https://github.com/inogs/OGSParaviewSuite/blob/master/OGSPlugins/OGSSpatialStats/doc/eq2.png" alt="" width="200"/>

This reduction of domains A and B requires to store the partial computations of the weights, mean and variance for each sub-domain:

<img src="https://github.com/inogs/OGSParaviewSuite/blob/master/OGSPlugins/OGSSpatialStats/doc/eq3.png" alt="" width="200"/>

A similar reduction can be implemented for the standard deviation:

<img src="https://github.com/inogs/OGSParaviewSuite/blob/master/OGSPlugins/OGSSpatialStats/doc/eq4.png" alt="" width="200"/>

the total standard derviation need to be divided by the sum of all the weights.

For the percentiles, it is necessary to order the values first (and retain the indices). Then, a percentile weight can be computed as:

<img src="https://github.com/inogs/OGSParaviewSuite/blob/master/OGSPlugins/OGSSpatialStats/doc/eq5.png" alt="" width="200"/>

where _'_ denotes the sorted weight. Then, the index of the value that is equal or immediately lower than the percentile weight must be obtained (e.g., for 50% the search is for wp = 0.5 or the imediately before). This index is called _s_, then (e.g., for 50%):

<img src="https://github.com/inogs/OGSParaviewSuite/blob/master/OGSPlugins/OGSSpatialStats/doc/eq6.png" alt="" width="200"/>