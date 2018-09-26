# Select basin, ParaView programable filter script.
#
# Programable filter intended to select a basin from 
# OGS data sets using the BasinMask variable.
#
# This filter is to be applied after applying the 
# "Select Coast" filter
#
# Basins are:
#	1:  Alboran Sea
#	2:  South Western Mediterranean (west)
#	3:  South Western Mediterranean (east)
#	4:  North Western Mediterranean
#	5:  Northern Tyrrhenian
#	6:  Southern Tyrrhenian
#	7:  Northern Adriatic
#	8:  Southern Adriatic
#	9:  Aegean Sea
#	10: Western Ionian
#	11: Eastern Ionian
#	12: Northern Ionian
#	13: Western Levantine
#	14: Northern Levantine
#	15: Southern Levantine
#	16: Eastern Levantine
#
# Combined form all the dataset.
#
# Arnau Miro, SoHPC 2017

Name  = 'OGSSelectBasin'
Label = 'OGS Select Basin'
Help  = ''

NumberOfInputs = 1
InputDataType  = 'vtkUnstructuredGrid'
OutputDataType = 'vtkUnstructuredGrid'
ExtraXml = ''

Properties = dict(
	Alboran_Sea                      = True, # alb  => 1
	South_Western_Mediterranean_west = True, # swm1 => 2
	South_Western_Mediterranean_east = True, # swm2 => 3
	North_Western_Mediterranean      = True, # nwm  => 4
	Northern_Tyrrhenian              = True, # tyr1 => 5
	Southern_Tyrrhenian              = True, # tyr2 => 6
	Northern_Adriatic                = True, # adr1 => 7
	Southern_Adriatic                = True, # adr2 => 8
	Aegean_Sea                       = True, # aeg  => 9
	Western_Ionian                   = True, # ion1 => 10
	Eastern_Ionian                   = True, # ion2 => 11
	Northern_Ionian                  = True, # ion3 => 12
	Western_Levantine                = True, # lev1 => 13
	Northern_Levantine               = True, # lev2 => 14
	Southern_Levantine               = True, # lev3 => 15
	Eastern_Levantine                = True  # lev4 => 16
	);

def RequestData():
	import vtk
	import numpy as np
	from vtk.util import numpy_support as npvtk

	# Get input data
	pdin = self.GetInput();

	# Recover the numpy array BasinMask from input
	BasinMask = npvtk.vtk_to_numpy( pdin.GetCellData().GetArray("basins mask") );

	# Create CutMask array full of zeros
	CutMask = np.zeros( np.shape(BasinMask) );

	# Fill with ones according to user input
	if (Alboran_Sea): CutMask[BasinMask == 1] = 1;
	if (South_Western_Mediterranean_west): CutMask[BasinMask == 2] = 1;
	if (South_Western_Mediterranean_east): CutMask[BasinMask == 3] = 1;
	if (North_Western_Mediterranean): CutMask[BasinMask == 4] = 1;
	if (Northern_Tyrrhenian): CutMask[BasinMask == 5] = 1;
	if (Southern_Tyrrhenian): CutMask[BasinMask == 6] = 1;
	if (Western_Ionian): CutMask[BasinMask == 10] = 1;
	if (Eastern_Ionian): CutMask[BasinMask == 11] = 1;
	if (Northern_Ionian): CutMask[BasinMask == 12] = 1;
	if (Northern_Adriatic): CutMask[BasinMask == 7] = 1;
	if (Southern_Adriatic): CutMask[BasinMask == 8] = 1;
	if (Western_Levantine): CutMask[BasinMask == 13] = 1;
	if (Northern_Levantine): CutMask[BasinMask == 14] = 1;
	if (Southern_Levantine): CutMask[BasinMask == 15] = 1;
	if (Eastern_Levantine): CutMask[BasinMask == 16] = 1;
	if (Aegean_Sea): CutMask[BasinMask == 9] = 1;

	# Add array to input
	vtkCutMask = npvtk.numpy_to_vtk(CutMask.ravel(),True,vtk.VTK_FLOAT);
	vtkCutMask.SetName("CutMask");
	pdin.GetCellData().AddArray(vtkCutMask);

	# Threshold construct
	thresh1 = vtk.vtkThreshold();
	thresh1.SetInputData(pdin);
	thresh1.ThresholdByUpper(1); # Erase whatever has 0
	thresh1.SetInputArrayToProcess(0, 0, 0, vtk.vtkDataObject.FIELD_ASSOCIATION_CELLS, "CutMask");
	thresh1.Update();

	# Grab the field and generate output
	field = thresh1.GetOutput();
	field.GetCellData().RemoveArray("CutMask"); # Remove the CutMask

	# Update the output port
	pdout = self.GetOutput().ShallowCopy(field);
