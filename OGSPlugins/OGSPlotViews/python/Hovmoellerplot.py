# Dummy file for the PythonView
#
# Arnau Miro, SoHPC 2017

Name  = 'Dummy'
Label = 'Dummy'
Help  = ''

NumberOfInputs = 1
InputDataType  = 'vtkUnstructuredGrid'
OutputDataType = 'vtkUnstructuredGrid'
ExtraXml = ''

Properties = dict(
	);

def RequestData():
	def setup_data(view):
		'''
		This function sets up the data to be plotted.
		'''
		view.EnableAllAttributeArrays()
	def render(view,width,height):
		'''
		This function decides the type of plot according 
		to the view type.

		Needs to be updated every time a plot view is created.
		'''
		view_type = view.GetClassName()
		# Vertical plot
		if view_type == "vtkOGSVerticalProfilePlot":
			return render_verticalplot(view, width, height)
		# Hovmoeller plot
		if view_type == "vtkOGSHovmoellerPlot":
			return render_hovmoellerplot(view, width, height)
		return None
	def smooth(a,WSZ):
		'''
		Smoothes a 1-D numpy array.
		
		WSZ: smoothing window size needs, which must be odd number,
		as in the original MATLAB implementation.
		'''
		import numpy as np
		out0  = np.convolve(a,np.ones(WSZ,dtype=int),'valid')/WSZ    
		r     = np.arange(1,WSZ-1,2)
		start = np.cumsum(a[:WSZ-1])[::2]/r
		stop  = (np.cumsum(a[:-WSZ:-1])[::2]/r)[::-1]
		return np.concatenate((  start , out0, stop  ))
	def render_hovmoellerplot(view, width, height):
		'''
		Performs Hovmoeller plots from data.
		'''
		from paraview import python_view
		from matplotlib import pyplot as plt
		
		import vtk, numpy as np
		from vtk.util import numpy_support as npvtk
		from datetime import datetime as dt

		# Create matplotlib figure
		figure = python_view.matplotlib_figure(width, height)

		# Set up the plot view
		ax = figure.add_subplot(1,1,1)
		# Title properties
		title_loc = 'center'
		if plot_title_alig == 0:
			title_loc = 'left'
		if plot_title_alig == 2:
			title_loc = 'right'
		ax.set_title(plot_title,
			fontsize=plot_title_font,
			fontweight='bold'    if plot_title_bold else None,
			style='italic'  if plot_title_ital else None,
			loc=title_loc)
		# X Axes properties
		ax.set_xlabel(x_label,
			fontsize=x_font,
			fontweight='bold' if x_bold else None,
			style='italic' if x_ital else None)
		# Y Axes properties
		ax.set_ylabel(y_label,
			fontsize=y_font,
			fontweight='bold' if y_bold else None,
			style='italic' if y_ital else None)
		if y_customrange: ax.set_ylim((y_min,y_max))

		# Plot data
		# Only get the first visible object, which should be a vtkTable
		if not view.GetNumberOfVisibleDataObjects() > 0: return None

		obj = view.GetVisibleDataObjectForRendering(0)
		if not obj: return None
		if (obj.GetClassName() != "vtkTable"): return None

		# Here we have mande sure we are dealing with a vtkTable
		# Now parse the columns and build the matrix data
		col_stack = []
		xIds      = []
		xNames    = []

		for cId in range(obj.GetNumberOfColumns()):
			# Recover the column
			column = obj.GetColumn(cId)
			# Deal with the depth column
			if column.GetName() == "depth":
				depth = -npvtk.vtk_to_numpy(column)
				continue
			# Else we have a data column
			xIds.append( cId )
			# Use datetime to format the date
			date = dt.strptime(column.GetName(),"%Y%m%d-%H:%M:%S")
			xNames.append( date.strftime(x_tick_form) )
			# Filter zeros and stack the column
			var = npvtk.vtk_to_numpy(column)
			var[var == 0.] = None
			# Algorithm to smooth data
			if smoothdata: var = smooth(var,smoothorder)
			col_stack.append( var )

		# Colormap selection
		cm = plt.cm.coolwarm
		if sel_colormap == 0: cm = plt.cm.jet
		if sel_colormap == 1: cm = plt.cm.rainbow
		if sel_colormap == 2: cm = plt.cm.hsv
		if sel_colormap == 3: cm = plt.cm.ocean
		if sel_colormap == 4: cm = plt.cm.terrain
		if sel_colormap == 5: cm = plt.cm.coolwarm
		if sel_colormap == 6: cm = plt.cm.viridis
		if sel_colormap == 7: cm = plt.cm.plasma
		if sel_colormap == 8: cm = plt.cm.inferno

		# Plot and colorbar
		ctf = ax.contourf(xIds,depth,np.column_stack(col_stack),cmap=cm)
		if draw_colorbar: figure.colorbar(ctf)

		# Set x axis
		ax.set_xticks(xIds)
		ax.set_xticklabels(xNames,rotation=x_rot)

		# Invert the Y axis
		ax.invert_yaxis()

		# Fit figure
		figure.tight_layout()

		# Save figure
		if savefigure and filename:
			figure.savefig(filename,dpi=outdpi,bbox_inches='tight')

		return python_view.figure_to_image(figure)