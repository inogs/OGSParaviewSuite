## OGS Omega Criterion

The **OGS Compute Omega Criterion** filter computes the Omega turbulence criterion from a _vtkRectilinearGrid_. The filter uses the same approach as the **OGS Derivatives** for the computation of the gradient.

The Omega vortex identification method (Liu, 2016) is a recent methjod proposed in order to overcome the scale difference of other criteria such as the Q-criterion or the lambda2 criterion. The Omega parameter is non-dimensional and fitted between 0 and 1. It is defined as:

<img src="https://github.com/inogs/OGSParaviewSuite/blob/master/OGSPlugins/OGSOmegaCriterion/doc/eq1.png" alt="" width="120"/>

where

<img src="https://github.com/inogs/OGSParaviewSuite/blob/master/OGSPlugins/OGSOmegaCriterion/doc/eq2.png" alt="" width="220"/>

and

<img src="https://github.com/inogs/OGSParaviewSuite/blob/master/OGSPlugins/OGSOmegaCriterion/doc/eq3.png" alt="" width="160"/>

The parameter _epsilon_ is a small enough parameter in order to avoid dividing by zero. Practically, it needs to be escaled with the problem in order to deal with the dimensionality. A correlation has been recently proposed (Dong, 2018) for _epsilon_ to adapt to the dimensions of the problem:

<img src="https://github.com/inogs/OGSParaviewSuite/blob/master/OGSPlugins/OGSOmegaCriterion/doc/eq4.png" alt="" width="180"/>

In the Mediterranean. a value of Omega = 0.5 captures the vortical structures in a similar manner than the region of Q < -Q0.

<img src="https://github.com/inogs/OGSParaviewSuite/blob/master/OGSPlugins/OGSOmegaCriterion/doc/Omega.png" alt="" width="1000"/>