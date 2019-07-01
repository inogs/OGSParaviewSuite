## OGS Derivatives

The **OGS Derivatives** filter is able to compute the spatial derivatives, i.e., the gradient of a signle scalar or vector array. This filter must be applied over a _RectilinearGrid_, as it uses a structured mesh approach for the gradient. Moreover, the filter can also compute the divergence, curl and Q of vectorial arrays.

The user can choose between five different gradient approaches:
* _FV 2nd Order_: Cell centered 2nd order finite volume approach.
* _FV 4th Order_: Cell centered 4th order finite volume approach.
* _OGSTM-BFM_: 1st order Euler approach as described in the NEMO book.
* _OGSTM-BFM2_: Experimental 2nd order central difference approach consistent with the NEMO book.
* _OGSTM-BFM4_: Experimental 4th order central difference approach consistent with the NEMO book.

By default, the _OGSTM-BFM_ method is used since it is consistent with the simulation.

### Gradient using Finite Volumes approach

Let a scalar field _q_ and a vector

<img src="https://github.com/inogs/OGSParaviewSuite/blob/master/OGSPlugins/OGSDerivatives/doc/eq1.png" alt="" width="160"/>

that live on the **T grid**. Their representation, using FV, is:

<img src="https://github.com/inogs/OGSParaviewSuite/blob/master/OGSPlugins/OGSDerivatives/doc/eq2.png" alt="" width="200"/>

The second order approach to the first derivative reads:

<img src="https://github.com/inogs/OGSParaviewSuite/blob/master/OGSPlugins/OGSDerivatives/doc/eq3.png" alt="" width="400"/>

for a scalar array or a component of a vectorial array. The derivations for the _y_ and _z_ components are analogous to these of the _x_ component.

<img src="https://github.com/inogs/OGSParaviewSuite/blob/master/OGSPlugins/OGSDerivatives/doc/eq4.png" alt="" width="400"/>

As mentioned, these derivatives live on the **T grid**, therefore, computing the divergence and curl is straightforward

<img src="https://github.com/inogs/OGSParaviewSuite/blob/master/OGSPlugins/OGSDerivatives/doc/eq5.png" alt="" width="400"/>

### Gradient using the NEMO approach

This implementation uses the approach described by [NEMO](https://www.nemo-ocean.eu/) and requires the mesh stretching tensors to be loaded by the **OGS Reader**.

Let a scalar field _q_ and a vector

<img src="https://github.com/inogs/OGSParaviewSuite/blob/master/OGSPlugins/OGSDerivatives/doc/eq1.png" alt="" width="160"/>

that live on the **T grid**. Their representation, according to NEMO documentation, is:

<img src="https://github.com/inogs/OGSParaviewSuite/blob/master/OGSPlugins/OGSDerivatives/doc/eq6.png" alt="" width="200"/>

The divergence of a vectorial field is:

<img src="https://github.com/inogs/OGSParaviewSuite/blob/master/OGSPlugins/OGSDerivatives/doc/eq7.png" alt="" width="200"/>

assuming that the derivatives of e1 and e2 are assumed to be very small or constant:

<img src="https://github.com/inogs/OGSParaviewSuite/blob/master/OGSPlugins/OGSDerivatives/doc/eq8.png" alt="" width="200"/>

Using discrete operators, the gradients of the scalar field and vectorial field are:

<img src="https://github.com/inogs/OGSParaviewSuite/blob/master/OGSPlugins/OGSDerivatives/doc/eq9.png" alt="" width="200"/>

where the delta operator basically denotes a forward Euler operation. Moreover, these gradients naturally fall in the **UVW grid**. The gradient can easily be expressed as:

<img src="https://github.com/inogs/OGSParaviewSuite/blob/master/OGSPlugins/OGSDerivatives/doc/eq10.png" alt="" width="200"/>

For the other boundary, the backwards difference is used. Two experimental gradient approaches are proposed using the second and fourth order central difference scheme previously illustrated. The divergence and the curl of a vectorial field that lives on the **T grid** is:

<img src="https://github.com/inogs/OGSParaviewSuite/blob/master/OGSPlugins/OGSDerivatives/doc/eq11.png" alt="" width="200"/>

Note that each of the components belong to the gradient living on the **T grid**. Therefore, the curl and the divergence are computed using the components of the gradient after this has been projected from the **UVW grid** to the **T grid**.
