## OGS Q-Criterion

The **OGS Compute Q-Criterion** computes the Q-criterion from a _vtkRectilinearGrid_ and generates a mask in a similar manner than the 
**OGS Compute Okubo-Weiss**. The filter uses the same approach as the **OGS Derivatives** for the computation of the gradient.

The Q-criterion (Hunt, 1989; Jeong and Hussain, 1995) is the three dimensional counterpart of the Okubo-Weiss parameter and is defined as:

<img src="https://github.com/inogs/OGSParaviewSuite/blob/master/OGSPlugins/OGSQCriterion/doc/eq1.png" alt="" width="150"/>

where:

<img src="https://github.com/inogs/OGSParaviewSuite/blob/master/OGSPlugins/OGSQCriterion/doc/eq2.png" alt="" width="300"/>

The relationship between the Okubo-Weiss parameter and the Q-criterion can be found, after some math, to be:

<img src="https://github.com/inogs/OGSParaviewSuite/blob/master/OGSPlugins/OGSQCriterion/doc/eq3.png" alt="" width="80"/>

This "correction" is applied so that the same three regions of the Okubo-Weiss parameter can be defined as:
* Vorticity dominated (Q < Q0): 0
* Background field (|Q| <= Q0): 1
* Strain dominated (Q > Q0): 2

with

<img src="https://github.com/inogs/OGSParaviewSuite/blob/master/OGSPlugins/OGSQCriterion/doc/eq4.png" alt="" width="80"/>

and 0.2 being a scaling factor.

This mask has the same properties as the Okubo-Weiss mask and can be managed by the **OGS Select Okubo-Weiss** filter by just changing the input mask field. As it can be seen, the Okubo-Weiss and the Q-criterion produce the same results for the surface of the Mediterranean sea.

<img src="https://github.com/inogs/OGSParaviewSuite/blob/master/OGSPlugins/OGSQCriterion/doc/QvsOW2.png" alt="" width="1000"/>