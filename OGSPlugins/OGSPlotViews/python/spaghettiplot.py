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
	def render_spaghettiplot(view, width, height):
		'''
		Performs spaghetti plots from data.
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
		if splot_title_alig == 0:
			title_loc = 'left'
		if splot_title_alig == 2:
			title_loc = 'right'
		if plt_vers > '1.3.1':
			ax.set_title(splot_title,
				fontsize=splot_title_font,
				fontweight='bold'    if splot_title_bold else None,
				style='italic'  if splot_title_ital else None,
				loc=title_loc)
		else:
			ax.set_title(splot_title,
				fontsize=splot_title_font,
				fontweight='bold'    if splot_title_bold else None,
				style='italic'  if splot_title_ital else None)			
		# X Axes properties
		ax.set_xlabel(sx_label,
			fontsize=sx_font,
			fontweight='bold' if sx_bold else None,
			style='italic' if sx_ital else None)
		# Y Axes properties
		ax.set_ylabel(sy_label,
			fontsize=sy_font,
			fontweight='bold' if sy_bold else None,
			style='italic' if sy_ital else None)
		if sy_customrange: ax.set_ylim((sy_min,sy_max))

		# Generate variable lists
		objlist = [s.split(';')[0] for s in labels.split('\n')] if not labels == '' else []

		xIdsg   = np.array([])
		xNamesg = np.array([])

		# Plot data
		if not view.GetNumberOfVisibleDataObjects() > 0: return python_view.figure_to_image(figure)
		for objId in range(view.GetNumberOfVisibleDataObjects()):
			# Get object that must be a vtkTable
			obj = view.GetVisibleDataObjectForRendering(objId)
			if not obj: continue
			if (obj.GetClassName() != "vtkTable"): continue
			# Here we have mande sure we are dealing with a vtkTable
			# Now parse the columns and build the array data
			var    = []
			xIds   = []
			xNames = []
			# Decide whether to plot the variable
			do_plot = False if len(objlist) > 0 else True
			args    = None
			kwargs  = None
			for ii in range(0,len(objlist)):
				if "%d" % (objId+1) == "%s" % (objlist[ii]): 
					do_plot = True
					if (len(labels.split('\n')[ii].split(';'))>1):
						args = labels.split('\n')[ii].split(';')[1]
					if (len(labels.split('\n')[ii].split(';'))>2):
						kwargs = dict(s.split('=') for s in labels.split('\n')[ii].split(';')[2].split(','))
			# Build the array data
			if do_plot:
				for cId in range(obj.GetNumberOfColumns()):
					# Recover the column
					column = obj.GetColumn(cId)
					# Else we have a data column
					xIds.append( cId )
					# Use datetime to format the date
					date = dt.strptime(column.GetName(),"%Y%m%d-%H:%M:%S")
					xNames.append( date.strftime(sx_tick_form) )
					# Filter zeros and stack the column
					val = npvtk.vtk_to_numpy(column)[0]
					val = None if val == 0. else val
					var.append(val)
				# Convert to numpy array
				z      = np.array(var)
				xIds   = np.array(xIds)
				xNames = np.array(xNames)
				# Smooth if requested
				if ssmoothdata: z = smooth(z,ssmoothorder)
				# Plot
				if args == None:
					ax.plot(xIds,z)
				elif kwargs == None:
					ax.plot(xIds,z,args)
				else:
					ax.plot(xIds,z,args,**kwargs)
				# Update global variables
				xIdsg   = xIds   if xIds.shape[0]   > xIdsg.shape[0]   else xIdsg
				xNamesg = xNames if xNames.shape[0] > xNamesg.shape[0] else xNamesg

		# Set x axis
		x_ticks = xIdsg if xIdsg.shape[0] < 10 else np.linspace(xIdsg[0],xIdsg[-1],sx_nticks,dtype=np.int)
		ax.set_xticks(x_ticks)
		ax.set_xticklabels(xNamesg[x_ticks],rotation=sx_rot)

		# Set legend
		if sshow_legend: 
			legend_prop = dict(size=slegend_font)
			if slegend_bold: legend_prop['weight'] = 'bold'
			if slegend_ital: legend_prop['style']  = 'italic'
			ax.legend(loc=slegend_location,prop=legend_prop)

		# Fit figure
		figure.tight_layout()

		# Save figure
		if ssavefigure and sfilename:
			# Save figure
			figure.savefig(sfilename,dpi=soutdpi,bbox_inches='tight')

		return python_view.figure_to_image(figure)