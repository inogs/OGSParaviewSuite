# Dummy file for the PythonView
#
# Arnau Miro, SoHPC 2017

# Variables


savefigure      = False
outname         = ''
outdpi          = 300

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
		Plot NETCDF using CARTOPY
		'''
		from paraview import python_view

		import MapPlotter as mp
		import matplotlib.pyplot as plt

		# Instance MapPlotter class
		plotter = mp.MapPlotter(projection='PlateCarree')

		# Set plot style
		if style == 0: plt.style.use('seaborn-white')
		if style == 1: plt.style.use('dark_background')
		if style == 2: plt.style.use('ggplot')

		# Create a figure
		figure = python_view.matplotlib_figure(width, height)
		ax     = figure.add_axes([0.08, 0.13, 0.78, 0.78,],projection=plotter.projection)

		# Title properties
		title_loc = 'center'
		if title_alig == 0: title_loc = 'left'
		if title_alig == 2: title_loc = 'right'

		# Set plot parameters
		params  = plotter.defaultParams()
		params['fig']       = figure
		params['ax']        = ax
		params['xlim']      = [x_min, x_max]
		params['ylim']      = [y_min, y_max]
		params['max_div']   = max_div
		params['grdstyle']  = {'size': grd_font, 'weight': 'bold' if grd_bold else None, 'style': 'italic' if grd_ital else None}
		params['features']  = []
		if show_coastline:  params['features'].append('coastline')
		if show_continents: params['features'].append('continents')
		if show_rivers:     params['features'].append('rivers')
		if show_image:      params['features'].append('image')
		params['img']       = imgfile
		params['title']     = [title,{'fontsize':title_font,'fontweight':'bold' if title_bold else None,'style':'italic' if title_ital else None,'loc':title_loc}]
		params['xlabel']    = [x_label,{'fontsize':x_font,'fontweight':'bold' if x_bold else None,'style':'italic' if x_ital else None}]
		params['ylabel']    = [y_label,{'fontsize':y_font,'fontweight':'bold' if y_bold else None,'style':'italic' if y_ital else None}]
		params['cmap']      = 'coolwarm'
		if sel_colormap == 0: params['cmap'] = 'jet'
		if sel_colormap == 1: params['cmap'] = 'rainbow'
		if sel_colormap == 2: params['cmap'] = 'hsv'
		if sel_colormap == 3: params['cmap'] = 'ocean'
		if sel_colormap == 4: params['cmap'] = 'terrain'
		if sel_colormap == 5: params['cmap'] = 'coolwarm'
		if sel_colormap == 6: params['cmap'] = 'viridis'
		if sel_colormap == 7: params['cmap'] = 'plasma'
		if sel_colormap == 8: params['cmap'] = 'inferno'
		params['bounds']    = [cbar_min,cbar_max]
		params['tick_font'] = grd_font
		params['label']     = {'label':cbar_label,'size':cbar_font,'weight':'bold' if cbar_bold else None,'style':'italic' if cbar_ital else None}

		# Plot
		if filename == '' or varname == '':
			figure = plotter.plot_empty(params=params)
		elif maskfile == '':
			figure = plotter.plot_from_file(filename,varname,lonname,latname,iTime=iTime,iDepth=iDepth,params=params)
		else:
			figure = plotter.plot_from_file_and_mask(filename,varname,maskfile,iTime=iTime,iDepth=iDepth,
				masklon=lonname,masklat=latname,params=params)

		# Save figure
		if savefigure and not outname == '':
			plotter.save(outname,dpi=outdpi)

		return python_view.figure_to_image(figure)