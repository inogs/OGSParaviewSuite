## OGS Writer

The **OGS Writer** lets the user export data into multiple formats. These formats are integrated to work with _vtkRectilinearGrid_, _vtkUnstructuredGrid_, _vtkPolyData_ and _vtkTable_ in both cell and point data format. All writers are time aware (i.e., they potentially allow to save a time series) with the exception of the _vtkTable_ writer. The available output formats are:

* _.npz_: Numpy compatible array dictionary for _vtkRectilinearGrid_, _vtkUnstructuredGrid_, _vtkPolyData_ and _vtkTable_. It can be read in python using **numpy**:
```python
import numpy as np

data = np.load('<EXAMPLE_NPZ>.npz')

data.keys() # Check arrays stored in data
data['xyz'] # Return xyz array from data dictionary
```

* _.field_: Binary array storage for _vtkRectilinearGrid_, _vtkUnstructuredGrid_ and _vtkPolyData_. This file must be read as follows:

|  type   |     Amount    | Name  |    Description    |
|---------|---------------|-------|-------------------|
| int32   |    1          | ncols | number of columns |
| int32   |    1          | mshsz | mesh size         |
| float32 |    1          | time  | simulation time   |
| float32 | 3 * mshsz     | xyz   | mesh points       |
| float32 | ncols * mshsz | data  | data array        |

* _.nc_: NetCDF-4 compatible file for _vtkRectilinearGrid_, _vtkUnstructuredGrid_ and _vtkPolyData_.
