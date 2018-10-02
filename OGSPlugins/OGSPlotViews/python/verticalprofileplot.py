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
	plot_title      = '',
	plot_title_font = 18,
	plot_title_bold = 0,
	plot_title_ital = 0,
	plot_title_alig = 1,
	x_label         = '',
	x_font          = 12,
	x_bold          = 0,
	x_ital          = 0,
	x_logscale      = 0,
	x_customrange   = 0,
	x_min           = 0.,
	x_max           = 0.,
	y_label         = '',
	y_font          = 12,
	y_bold          = 0,
	y_ital          = 0,
	y_logscale      = 0,
	y_customrange   = 0,
	y_min           = 0.,
	y_max           = 0.,
	y_fact          = 1000.,
	minorticks      = 0,
	show_grid       = 0,
	show_legend     = 0,
	filterzeros     = 0,
	legend_location = 0,
	legend_font     = 0,
	legend_bold     = 0,
	legend_ital     = 0
	);

def RequestData():
	def setup_data(view):
		view.EnableAllAttributeArrays()
	def render(view, width, height):
		from paraview import python_view
		figure = python_view.matplotlib_figure(width, height)
		import vtk
		from vtk.util import numpy_support as npvtk

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
		if x_logscale: ax.set_xscale('log')
		if x_customrange: ax.set_xlim((x_min,x_max))
		# Y Axes properties
		ax.set_ylabel(y_label,
			fontsize=y_font,
			fontweight='bold' if y_bold else None,
			style='italic' if y_ital else None)
		if y_logscale: ax.set_yscale('log')
		if y_customrange: ax.set_ylim((y_min,y_max))

		if minorticks: ax.minorticks_on()
		if show_grid: ax.grid()

		# Generate variable lists
		objlist = [s.split(';')[1] for s in variables.split('\n')] if not variables == '' else []
		varlist = [s.split(';')[0] for s in variables.split('\n')] if not variables == '' else []

		# Plot data
		for objId in xrange(view.GetNumberOfVisibleDataObjects()):
			obj = view.GetVisibleDataObjectForRendering(objId)
			# Recover the depth
			z   = npvtk.vtk_to_numpy(obj.GetPoints().GetData())[:,2]/y_fact
			# For each variable, decide if the plot is needed
			for varId in xrange(view.GetNumberOfAttributeArrays(objId,vtk.vtkDataObject.POINT)):         
				varname = view.GetAttributeArrayName(objId,vtk.vtkDataObject.POINT,varId)
				# Filter some variables (Do not plot them)
				if varname == 'basins mask':       continue
				if varname == 'coast mask':        continue
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
					if filterzeros: var[var == 0.] = None
					if args == None:
						ax.plot(var,-z,label=varname)
					elif kwargs == None:
						ax.plot(var,-z,args,label=varname)
					else:
						ax.plot(var,-z,args,**kwargs)

		# Invert the Y axis
		ax.invert_yaxis()

		# Set legend
		if show_legend: 
			legend_prop = dict(size=legend_font)
			if legend_bold: legend_prop['weight'] = 'bold'
			if legend_ital: legend_prop['style']  = 'italic'
			ax.legend(loc=legend_location,prop=legend_prop)

		return python_view.figure_to_image(figure)