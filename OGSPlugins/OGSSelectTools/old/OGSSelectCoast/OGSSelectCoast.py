# Select coast, ParaView programable filter script.
#
# Programable filter intended to select a coast type from 
# OGS data sets using the CoastMask variable. Afterwards
# removes the CoastMask array since it is no longer needed.
#
# This filter is to be applied after loading the data or
# computing the OkuboWeiss criterion.
#
# Coast values are:
#	1: coast
#	2: open sea
#
# Combined form all the dataset.
#
# Arnau Miro, SoHPC 2017

Name  = 'OGSSelectCoast'
Label = 'OGS Select Coast'
Help  = ''

NumberOfInputs = 1
InputDataType  = 'vtkRectilinearGrid'
OutputDataType = 'vtkUnstructuredGrid'
ExtraXml = ''

Properties = dict(
	coast = True,    # Value of 1 in CoastMask
	open_sea = True  # Value of 2 in CoastMask
	);

def RequestData():
	import vtk
	import numpy as np
	from vtk.util import numpy_support as npvtk

	# Get input data
	pdin = self.GetInput();

	# Threshold construct
	thresh1 = vtk.vtkThreshold();
	thresh1.SetInputData(pdin);

	# Apply a threshold filter according to the CoastMask
	if coast and open_sea:     thresh1.ThresholdBetween(1,2);
	if coast and not open_sea: thresh1.ThresholdBetween(1,1.5);
	if not coast and open_sea: thresh1.ThresholdBetween(1.5,2);
	if not coast and not open_sea: thresh1.ThresholdBetween(0,0.5); # Just the ground

	# Apply threshold on CoastMask
	thresh1.SetInputArrayToProcess(0, 0, 0, vtk.vtkDataObject.FIELD_ASSOCIATION_CELLS, "coast mask");
	# Update filter
	thresh1.Update();

	# Grab the vtkField from threshold filter
	field = thresh1.GetOutput();
	field.GetCellData().RemoveArray("coast mask"); # Remove the coast mask (no longer needed)

	# Update the output port
	pdout = self.GetOutput().ShallowCopy(field);