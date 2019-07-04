## Utilities

These are a collection of utility and helper codes that complement the **OGS ParaView Suite**. They are thoroughly used throughout the filters. Some of these tools are described here.

### Field and VTK interface

The template class field::Field is an optimized container for multidimensional mesh array variables (e.g., velocity). It contains an iterator and works in a much similar way to the std::vector class. It is declared as:

```cpp
field::Field<type> f(n,m)
```

where:

* _n_: number of dimensions on the mesh.
* _m_: number of dimensions of the array.

To recover a value use:

```cpp
f[idx_n][idx_m]
```

where:

* _idx_n_: is the index associated to _n_.
* _idx_m_: is the index associated to _m_.

More information can be found in the _field.h_ file regarding the methods and operators defined, as well as the iterator. Template methods to switch from _Field_ to _VTK_ are provided in _vtkFields.hpp_.

```cpp
VTKARRAY createVTKfromField<VTKARRAY,type>(name, field)
```

is an easy interface to generate _VTK_ arrays from _Field_ variables.

```cpp
field::Field<type> createFieldfromVTK<VTKARRAY,type>(vtkArray)
```

is an easy interface to generate _Field_ from _VTK_ arrays.

### Conversion from UVW grid to T grid

Suppose a generic vector

<img src="https://github.com/inogs/OGSParaviewSuite/blob/master/OGSPlugins/_utils/doc/eq1.png" alt="" width="150"/>

that falls on the **UVW grid**. In order to obtain its representation on the **T grid**, let us consider the fluxes on the i-1 and i faces of the cell (for the first component of the vector):

<img src="https://github.com/inogs/OGSParaviewSuite/blob/master/OGSPlugins/_utils/doc/eq2.png" alt="" width="220"/>

The flux in the middle of the cell, i.e., T is:

<img src="https://github.com/inogs/OGSParaviewSuite/blob/master/OGSPlugins/_utils/doc/eq3.png" alt="" width="150"/>

Let's assume that the flux in the middle of the cell can be obtained by averaging that of the faces:

<img src="https://github.com/inogs/OGSParaviewSuite/blob/master/OGSPlugins/_utils/doc/eq4.png" alt="" width="150"/>

Therefore:

<img src="https://github.com/inogs/OGSParaviewSuite/blob/master/OGSPlugins/_utils/doc/eq5.png" alt="" width="400"/>

In a similar manner, the other components can be obtained as:

<img src="https://github.com/inogs/OGSParaviewSuite/blob/master/OGSPlugins/_utils/doc/eq6.png" alt="" width="380"/>

Since, for the third component, e1w = e1t and e2w = e2t.

A code that implements the conversion from the **UVW grid** to the **T grid** is included in _fieldOperations.hpp_:

```cpp
field::Field<type> UVW2T(field_UVW, e1, e2, e3, nx, ny, nz)
```