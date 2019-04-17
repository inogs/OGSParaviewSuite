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
		# Spaghetti plot
		if view_type == "vtkOGSSpaghettiPlot":
			return render_spaghettiplot(view, width, height)
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
		from matplotlib import __version__ as plt_vers
		from matplotlib import pyplot as plt
		from matplotlib.ticker import MaxNLocator
		
		import vtk, numpy as np
		from vtk.util import numpy_support as npvtk
		from datetime import datetime as dt

		# Create matplotlib figure
		figure = python_view.matplotlib_figure(width, height)

		# Set up the plot view
		ax = figure.add_subplot(1,1,1)
		# Title properties
		title_loc = 'center'
		if hplot_title_alig == 0:
			title_loc = 'left'
		if hplot_title_alig == 2:
			title_loc = 'right'
		# Fix concerning matplotlib 1.1.1 that is shipped with the superbuild
		if plt_vers > '1.3.1':
			ax.set_title(hplot_title,
				fontsize=hplot_title_font,
				fontweight='bold'    if hplot_title_bold else None,
				style='italic'  if hplot_title_ital else None,
				loc=title_loc)
		else:
			ax.set_title(hplot_title,
				fontsize=hplot_title_font,
				fontweight='bold'    if hplot_title_bold else None,
				style='italic'  if hplot_title_ital else None)
		# X Axes properties
		ax.set_xlabel(hx_label,
			fontsize=hx_font,
			fontweight='bold' if hx_bold else None,
			style='italic' if hx_ital else None)
		# Y Axes properties
		ax.set_ylabel(hy_label,
			fontsize=hy_font,
			fontweight='bold' if hy_bold else None,
			style='italic' if hy_ital else None)
		if hy_customrange: ax.set_ylim((hy_min,hy_max))

		# Plot data
		# Only get the first visible object, which should be a vtkTable
		if not view.GetNumberOfVisibleDataObjects() > 0: return python_view.figure_to_image(figure)

		obj = view.GetVisibleDataObjectForRendering(0)
		if not obj: return python_view.figure_to_image(figure)
		if (obj.GetClassName() != "vtkTable"): return python_view.figure_to_image(figure)

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
			xNames.append( date.strftime(hx_tick_form) )
			# Filter zeros and stack the column
			var = npvtk.vtk_to_numpy(column)
			var[var == 0.] = None
			# Algorithm to smooth data
			if hsmoothdata: var = smooth(var,hsmoothorder)
			col_stack.append( var )

		z      = np.column_stack(col_stack)
		xIds   = np.array(xIds)
		xNames = np.array(xNames)

		# Colormap selection
		cm = plt.cm.coolwarm
		if hsel_colormap == 0: cm = plt.cm.jet
		if hsel_colormap == 1: cm = plt.cm.rainbow
		if hsel_colormap == 2: cm = plt.cm.hsv
		if hsel_colormap == 3: cm = plt.cm.ocean
		if hsel_colormap == 4: cm = plt.cm.terrain
		if hsel_colormap == 5: cm = plt.cm.coolwarm
		if hsel_colormap == 6: cm = plt.cm.viridis
		if hsel_colormap == 7: cm = plt.cm.plasma
		if hsel_colormap == 8: cm = plt.cm.inferno

		# Plot and colorbar
		z_min = np.nanmin(z)
		z_max = np.nanmax(z)
		ctf = ax.contourf(xIds,depth,z,cmap=cm,levels=np.linspace(z_min,z_max,256))
		if hdraw_colorbar: 
			cbar = figure.colorbar(ctf)
			cbar.set_ticks(np.linspace(z_min,z_max,10))

		# Set x axis
		x_ticks = xIds if xIds.shape[0] < 10 else np.linspace(xIds[0],xIds[-1],hx_nticks,dtype=np.int)
		ax.set_xticks(x_ticks)
		ax.set_xticklabels(xNames[x_ticks],rotation=hx_rot)

		# Invert the Y axis
		ax.invert_yaxis()

		# Fit figure
		figure.tight_layout()

		# Save figure
		if hsavefigure and hfilename:
			# Fix ugly white lines in PDF
			for c in ctf.collections: c.set_edgecolor("face")
			# Save figure
			figure.savefig(hfilename,dpi=houtdpi,bbox_inches='tight')

		return python_view.figure_to_image(figure)