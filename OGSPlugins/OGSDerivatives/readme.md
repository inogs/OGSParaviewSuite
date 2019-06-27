## OGS Derivatives

The **OGS Derivatives** filter is able to compute the spatial derivatives, i.e., the gradient of a signle scalar or vector array. This filter must be applied over a _RectilinearGrid_, as it uses a structured mesh approach for the gradient. Moreover, the filter can also compute the divergence, curl and Q of vectorial arrays.

The user can choose between five different gradient approaches:
* _FV 2nd Order_: Cell centered 2nd order finite volume approach.
* _FV 4th Order_: Cell centered 4th order finite volume approach.
* _OGSTM-BFM_: 1st order Euler approach as described in the NEMO book.
* _OGSTM-BFM2_: Experimental 2nd order central difference approach consistent with the NEMO book.
* _OGSTM-BFM4_: Experimental 4th order central difference approach consistent with the NEMO book.

By default, the _OGSTM-BFM_ method is used since it is consistent with the simulation.