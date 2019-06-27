## OGS Gradient

The **OGS Gradient** filter computes the spatial derivatives of a field from any kind of grid. It is based on the _ParaView_ filter _Gradient of Unstructured DataSet_, therefore it uses _ParaView_'s algorithm for the computation of the gradient. The filter is also able to return the vorticity (curl) and divergence of vectorial fields, as well as the Q, Lambda2, and Omega criteria for turbulence.

The gradient algorithm has been modified to be consistent with the scaling on the _z_ direction imposed in the **OGS Reader**. This is achieved by multiplying the gradient in the _z_ direction by the scaling factor. The gradients obtained with this filter are shown to be consistent with the ones obtained by the **OGS Derivatives**.