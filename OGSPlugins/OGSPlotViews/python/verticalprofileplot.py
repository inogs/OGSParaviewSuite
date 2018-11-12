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
	def render_verticalplot(view, width, height):
		'''
		Performs vertical plots from data.
		'''
		from paraview import python_view
		from matplotlib import pyplot as plt
		
		import vtk, numpy as np
		from vtk.util import numpy_support as npvtk

		# Create matplotlib figure
		figure = python_view.matplotlib_figure(width, height)

		# Use xkcd
		if pplot_xkcd: plt.xkcd() 
		else:         plt.rcdefaults()

		# Set up the plot view
		ax = figure.add_subplot(1,1,1)
		# Title properties
		title_loc = 'center'
		if pplot_title_alig == 0:
			title_loc = 'left'
		if pplot_title_alig == 2:
			title_loc = 'right'
		ax.set_title(pplot_title,
			fontsize=pplot_title_font,
			fontweight='bold'    if pplot_title_bold else None,
			style='italic'  if pplot_title_ital else None,
			loc=title_loc)
		# X Axes properties
		ax.set_xlabel(px_label,
			fontsize=px_font,
			fontweight='bold' if px_bold else None,
			style='italic' if px_ital else None)
		if px_logscale: ax.set_xscale('log')
		if px_customrange: ax.set_xlim((px_min,px_max))
		# Y Axes properties
		ax.set_ylabel(py_label,
			fontsize=py_font,
			fontweight='bold' if py_bold else None,
			style='italic' if py_ital else None)
		if py_logscale: ax.set_yscale('log')
		if py_customrange: ax.set_ylim((py_min,py_max))

		if pminorticks: ax.minorticks_on()
		if pshow_grid: ax.grid()

		# Generate variable lists
		objlist = [s.split(';')[1] for s in variables.split('\n')] if not variables == '' else []
		varlist = [s.split(';')[0] for s in variables.split('\n')] if not variables == '' else []

		# Plot data
		for objId in xrange(view.GetNumberOfVisibleDataObjects()):
			obj = view.GetVisibleDataObjectForRendering(objId)
			# Recover the depth
			z   = npvtk.vtk_to_numpy(obj.GetPoints().GetData())[:,2]/py_fact
			# For each variable, decide if the plot is needed
			for varId in xrange(view.GetNumberOfAttributeArrays(objId,vtk.vtkDataObject.POINT)):         
				varname = view.GetAttributeArrayName(objId,vtk.vtkDataObject.POINT,varId)
				# Filter some variables (Do not plot them)
				if varname == 'basins mask':       continue
				if varname == 'coast mask':        continue
				if varname == 'e1':                continue
				if varname == 'e2':                continue
				if varname == 'e3':                continue
				if varname == 'vtkValidPointMask': continue
				# Decide whether to plot the variable
				do_plot = False if len(varlist) > 0 else True
				args    = None
				kwargs  = None
				for ii in xrange(0,len(varlist)):
					if "%s %d" % (varname,objId+1) == "%s %s" % (varlist[ii],objlist[ii]): 
						do_plot = True
						if (len(variables.split('\n')[ii].split(';'))>2):
							args = variables.split('\n')[ii].split(';')[2]
						if (len(variables.split('\n')[ii].split(';'))>3):
							kwargs = dict(s.split('=') for s in variables.split('\n')[ii].split(';')[3].split(','))
				# Plot
				if do_plot:
					var = npvtk.vtk_to_numpy(obj.GetPointData().GetArray(varname))
					# Filter zeros inside 
					if pfilterzeros: var[var == 0.] = None
					# Algorithm to smooth data
					if psmoothdata: var = smooth(var,psmoothorder)
					# Plot
					if args == None:
						ax.plot(var,-z,label=varname)
					elif kwargs == None:
						ax.plot(var,-z,args,label=varname)
					else:
						ax.plot(var,-z,args,**kwargs)

		# Invert the Y axis
		ax.invert_yaxis()

		# Set legend
		if pshow_legend: 
			legend_prop = dict(size=plegend_font)
			if plegend_bold: legend_prop['weight'] = 'bold'
			if plegend_ital: legend_prop['style']  = 'italic'
			ax.legend(loc=plegend_location,prop=legend_prop)

		if psavefigure and pfilename:
			figure.savefig(pfilename,dpi=poutdpi,bbox_inches='tight')

		return python_view.figure_to_image(figure)