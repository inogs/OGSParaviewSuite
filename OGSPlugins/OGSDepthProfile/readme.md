## OGS Depth Profile

The **OGS Depth Profile** interpolates a line on the mesh and generates a profile that is stored as a *vtkpolyline*, to be plotted on the **OGS Vertical Profile Plot** of thje **OGS Plot Tools**. 

This plugin can either use a high resolution line source or a **OGS Depth Line Source** from the **OGS Utils** as input and the dataset as source; which can be either cell or point data.
* For _CellData_, the interpolation takes the closest calue of the cell to the line point. No further interpolation is performed.
* For _PointData_, it uses the interpolate tuple method to interpolate the value from the points of the closest cell using a weighted interpolation.

Arrays such as the mask arrays or the mesh stretching factors are ignored when computing the profiles.