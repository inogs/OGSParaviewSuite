# OGS Vortex Identification

This plugin provides a series of tools to identify vortices. A list of each one of them is provided below:
* Okubo-Weiss parameter.
* Q-criterion.
* Lambda2 criterion.
* Omega criterion.

## OGS Okubo-Weiss

The **OGS Compute Okubo-Weiss** and the **OGS Select Okubo-Weiss** are two filters that deal with the computation of the Okubo-Weiss parameter and eddy identification. The filter uses the same approach as the **OGS Derivatives** for the computation of the gradient.

The Okubo-Weiss parameter (Okubo, 1970; Weiss, 1991) is computed as:

<img src="https://github.com/inogs/OGSParaviewSuite/blob/master/OGSPlugins/OGSVortexIdentification/doc/OWeq1.png" alt="" width="160"/>

where:

<img src="https://github.com/inogs/OGSParaviewSuite/blob/master/OGSPlugins/OGSVortexIdentification/doc/OWeq2.png" alt="" width="250"/>

An eddy is then localized as a vorticity dominated area (or eddy core) surrounded by a strain dominated region (circulation cell) as defined by (Isern-Fontanet. 2004). Then:

<img src="https://github.com/inogs/OGSParaviewSuite/blob/master/OGSPlugins/OGSVortexIdentification/doc/OW.png" alt="" width="500"/>

with:

<img src="https://github.com/inogs/OGSParaviewSuite/blob/master/OGSPlugins/OGSVortexIdentification/doc/OWeq3.png" alt="" width="120"/>

and 0.2 being a scaling factor.

The **OGS Compute Okubo-Weiss** also computes a mask that can be used to select different regions using the **OGS Select Okubo-Weiss**. This mask is constructed so that:
* Vorticity dominated: 0
* Background field: 1
* Strain dominated: 2

This convention is used instead of -1, 0 and 1 so that the array can be stored as an _uint8_.

## OGS Q-Criterion

The **OGS Compute Q-Criterion** computes the Q-criterion from a _vtkRectilinearGrid_ and generates a mask in a similar manner than the 
**OGS Compute Okubo-Weiss**. The filter uses the same approach as the **OGS Derivatives** for the computation of the gradient.

The Q-criterion (Hunt, 1989; Jeong and Hussain, 1995) is the three dimensional counterpart of the Okubo-Weiss parameter and is defined as:

<img src="https://github.com/inogs/OGSParaviewSuite/blob/master/OGSPlugins/OGSVortexIdentification/doc/Qeq1.png" alt="" width="150"/>

where:

<img src="https://github.com/inogs/OGSParaviewSuite/blob/master/OGSPlugins/OGSVortexIdentification/doc/Qeq2.png" alt="" width="300"/>

The relationship between the Okubo-Weiss parameter and the Q-criterion can be found, after some math, to be:

<img src="https://github.com/inogs/OGSParaviewSuite/blob/master/OGSPlugins/OGSVortexIdentification/doc/Qeq3.png" alt="" width="80"/>

This "correction" is applied so that the same three regions of the Okubo-Weiss parameter can be defined as:
* Vorticity dominated (Q < Q0): 0
* Background field (|Q| <= Q0): 1
* Strain dominated (Q > Q0): 2

with

<img src="https://github.com/inogs/OGSParaviewSuite/blob/master/OGSPlugins/OGSVortexIdentification/doc/Qeq4.png" alt="" width="80"/>

and 0.2 being a scaling factor.

This mask has the same properties as the Okubo-Weiss mask and can be managed by the **OGS Select Okubo-Weiss** filter by just changing the input mask field. As it can be seen, the Okubo-Weiss and the Q-criterion produce the same results for the surface of the Mediterranean sea.

<img src="https://github.com/inogs/OGSParaviewSuite/blob/master/OGSPlugins/OGSVortexIdentification/doc/QvsOW2.png" alt="" width="1000"/>

## OGS Lambda2 Criterion

The **OGS Compute Lambda2 Criterion** computes the lambda2 turbulence criterion from a _vtkRectilinearGrid_ (structured grid) and generates the lambda2 value. The filter uses the same approach as the **OGS Derivatives** for the computation of the gradient.

The lambda2 criterion (Jeong & Hussain, 1995) is based on computing the value of the second eigenvalue (lambda2) of the tensor S2 + O2 where:

<center><img src="https://github.com/inogs/OGSParaviewSuite/blob/master/OGSPlugins/OGSVortexIdentification/doc/L2eqs.png" alt="" width="200"/></center>

Properly tuned, this criterion yelds the same result as the Q criterion.

The code makes use of the _matrixMN class_ to deal with matrix operations and the computation of the eigenvalues. For the latter, a Jacobi iteration algorithm is used.

## OGS Omega Criterion

The **OGS Compute Omega Criterion** filter computes the Omega turbulence criterion from a _vtkRectilinearGrid_. The filter uses the same approach as the **OGS Derivatives** for the computation of the gradient.

The Omega vortex identification method (Liu, 2016) is a recent methjod proposed in order to overcome the scale difference of other criteria such as the Q-criterion or the lambda2 criterion. The Omega parameter is non-dimensional and fitted between 0 and 1. It is defined as:

<img src="https://github.com/inogs/OGSParaviewSuite/blob/master/OGSPlugins/OGSVortexIdentification/doc/Oeq1.png" alt="" width="120"/>

where

<img src="https://github.com/inogs/OGSParaviewSuite/blob/master/OGSPlugins/OGSVortexIdentification/doc/Oeq2.png" alt="" width="220"/>

and

<img src="https://github.com/inogs/OGSParaviewSuite/blob/master/OGSPlugins/OGSVortexIdentification/doc/Oeq3.png" alt="" width="160"/>

The parameter _epsilon_ is a small enough parameter in order to avoid dividing by zero. Practically, it needs to be escaled with the problem in order to deal with the dimensionality. A correlation has been recently proposed (Dong, 2018) for _epsilon_ to adapt to the dimensions of the problem:

<img src="https://github.com/inogs/OGSParaviewSuite/blob/master/OGSPlugins/OGSVortexIdentification/doc/Oeq4.png" alt="" width="180"/>

In the Mediterranean. a value of Omega > 0.5 captures the vortical structures in a similar manner than the region of Q < -Q0.

<img src="https://github.com/inogs/OGSParaviewSuite/blob/master/OGSPlugins/OGSVortexIdentification/doc/Omega.png" alt="" width="1000"/>