## OGS Okubo-Weiss

The **OGS Compute Okubo-Weiss** and the **OGS Select Okubo-Weiss** are two filters that deal with the computation of the Okubo-Weiss parameter and eddy identification. The filter uses the same approach as the **OGS Derivatives** for the computation of the gradient.

The Okubo-Weiss parameter (Okubo, 1970; Weiss, 1991) is computed as:

<img src="https://github.com/inogs/OGSParaviewSuite/blob/master/OGSPlugins/OGSOkuboWeiss/doc/eq1.png" alt="" width="160"/>

where:

<img src="https://github.com/inogs/OGSParaviewSuite/blob/master/OGSPlugins/OGSOkuboWeiss/doc/eq2.png" alt="" width="250"/>

An eddy is then localized as a vorticity dominated area (or eddy core) surrounded by a strain dominated region (circulation cell) as defined by (Isern-Fontanet. 2004). Then:

<img src="https://github.com/inogs/OGSParaviewSuite/blob/master/OGSPlugins/OGSOkuboWeiss/doc/OW.png" alt="" width="500"/>

with:

<img src="https://github.com/inogs/OGSParaviewSuite/blob/master/OGSPlugins/OGSOkuboWeiss/doc/eq3.png" alt="" width="120"/>

and 0.2 being a scaling factor.

The **OGS Compute Okubo-Weiss** also computes a mask that can be used to select different regions using the **OGS Select Okubo-Weiss**. This mask is constructed so that:
* Vorticity dominated: 0
* Background field: 1
* Strain dominated: 2

This convention is used instead of -1, 0 and 1 so that the array can be stored as an _uint8_.