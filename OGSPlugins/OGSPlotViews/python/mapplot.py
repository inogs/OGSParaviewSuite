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
		# MapPlot plot
		if view_type == "vtkOGSMapPlot":
			return render_mapplot(view, width, height)
		return None
	def render_mapplot(view,width,height):
		'''
		Performs map plots from slice data.
		'''
		from paraview import python_view

		import io, requests
		import cartopy.crs as ccrs, cartopy.feature as cfeature
		import vtk, numpy as np, matplotlib, matplotlib.pyplot as plt

		from vtk.util import numpy_support as npvtk
		from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

		# Select projection
		projection = ccrs.PlateCarree() # muse_projection
		
		# Set plot style
		if mstyle == 0: plt.style.use('seaborn-white')
		if mstyle == 1: plt.style.use('dark_background')
		if mstyle == 2: plt.style.use('ggplot')

		# Create figure and axes
		figure = python_view.matplotlib_figure(width, height)
		ax     = figure.add_axes([0.08, 0.13, 0.78, 0.78,],projection=projection)

		# Draw coastline
		if mshow_coastline:
			ax.add_feature(
					cfeature.NaturalEarthFeature('physical', 'coastline', '50m', 
						edgecolor=(0,0,0,1), facecolor='none', linewidth=1.)
				)
		if mshow_borders:
			ax.add_feature(
					cfeature.NaturalEarthFeature('cultural', 'admin_0_boundary_lines_land', '50m', 
						edgecolor=(0,0,0,.6), facecolor='none', linewidth=.95)
				)
		if mshow_rivers:
			ax.add_feature(
					cfeature.NaturalEarthFeature('physical', 'rivers_lake_centerlines', '50m', 
						edgecolor=(0.0,0.0,0.545,.75), facecolor='none', linewidth=.75)
				)

		if mshow_image and mimg_file == '':
			ax.stock_img()
		if mshow_image and not mimg_file == '':
			# Detect if we are dealing with a URL or a path
			if 'https://' in img or 'http://' in img:
				ax.imshow(plt.imread(io.BytesIO(requests.get(img).content)), origin='upper', 
					transform=ccrs.PlateCarree(), extent=[-180, 180, -90, 90])
			else:
				ax.imshow(plt.imread(img), origin='upper', 
					transform=ccrs.PlateCarree(), extent=[-180, 180, -90, 90])

		# Title properties
		title_loc = 'center'
		if mplot_title_alig == 0:
			title_loc = 'left'
		if mplot_title_alig == 2:
			title_loc = 'right'
		ax.set_title(mplot_title,
			fontsize=mplot_title_font,
			fontweight='bold'    if mplot_title_bold else None,
			style='italic'  if mplot_title_ital else None,
			loc=title_loc)

		# X Axes properties
		ax.set_xlabel(mx_label,
			fontsize=mx_font,
			fontweight='bold' if mx_bold else None,
			style='italic' if mx_ital else None)
		# Y Axes properties
		ax.set_ylabel(my_label,
			fontsize=my_font,
			fontweight='bold' if my_bold else None,
			style='italic' if my_ital else None)

		# Set axis limits
		xlim = [mx_min, mx_max]; ax.set_xlim(xlim)
		ylim = [my_min, my_max]; ax.set_ylim(ylim)
		
		gl = ax.gridlines(crs=projection,draw_labels=True,linewidth=0)
		gl.xlabels_top   = False
		gl.ylabels_right = False
		gl.xlocator      = matplotlib.ticker.FixedLocator(np.arange(xlim[0],xlim[1],max_div))
		gl.ylocator      = matplotlib.ticker.FixedLocator(np.arange(ylim[0],ylim[1],max_div))
		gl.xformatter    = LONGITUDE_FORMATTER
		gl.yformatter    = LATITUDE_FORMATTER
		gl.xlabel_style  = {'size': mx_font, 'weight': 'bold' if mx_bold else None, 'style': 'italic' if mx_ital else None}
		gl.ylabel_style  = {'size': my_font, 'weight': 'bold' if my_bold else None, 'style': 'italic' if my_ital else None}

		# Filter data type
		if not view.GetNumberOfVisibleDataObjects() > 0: return python_view.figure_to_image(figure)

		obj = view.GetVisibleDataObjectForRendering(0)
		if not obj: return python_view.figure_to_image(figure)
		
		if (obj.GetClassName() != "vtkPolyData"): return python_view.figure_to_image(figure)

		# Only work with point data
		if (obj.GetPointData().GetNumberOfArrays() == 0): return python_view.figure_to_image(figure)

		# Check loaded variable
		if not obj.GetPointData().GetArray(mvarname): return python_view.figure_to_image(figure)

		# Load data arrays
		metadata = obj.GetFieldData().GetAbstractArray('Metadata')
		xyz      = npvtk.vtk_to_numpy(obj.GetPoints().GetData())
		data     = npvtk.vtk_to_numpy(obj.GetPointData().GetArray(mvarname))

		# Deal with vectorial arrays
		if len(data.shape) == 2:
			if not msel_comp in [0,1,2]:
				data = np.linalg.norm(data,axis=1)
			else:
				data = data[:,msel_comp]

		# Obtain the triangulation from VTK
		triang = [[int(obj.GetCell(icell).GetPointId(p)) for p in range(obj.GetCell(icell).GetNumberOfPoints())] 
			for icell in range(obj.GetNumberOfCells())]

		# Colormap
		cmap = plt.get_cmap('coolwarm',256)
		if msel_colormap == 0: cmap = plt.get_cmap('jet',256)
		if msel_colormap == 1: cmap = plt.get_cmap('rainbow',256)
		if msel_colormap == 2: cmap = plt.get_cmap('hsv',256)
		if msel_colormap == 3: cmap = plt.get_cmap('ocean',256)
		if msel_colormap == 4: cmap = plt.get_cmap('terrain',256)
		if msel_colormap == 5: cmap = plt.get_cmap('coolwarm',256)
		if msel_colormap == 6: cmap = plt.get_cmap('viridis',256)
		if msel_colormap == 7: cmap = plt.get_cmap('plasma',256)
		if msel_colormap == 8: cmap = plt.get_cmap('inferno',256)

		# XYZ must be projected to the projection used
		pv_projection = ccrs.PlateCarree(); 
		proj_name = metadata.GetValue(7)
		if proj_name == 'Mercator':     pv_projection = ccrs.Mercator(false_easting=989634.3811336625, false_northing=-3512473.95569)
		if proj_name == 'Cylindrical':  pv_projection = ccrs.epsg(4087)
		if proj_name == 'Google':       pv_projection = ccrs.GOOGLE_MERCATOR
		if proj_name == 'Mollweide':    pv_projection = ccrs.Mollweide()
		if proj_name == 'Orthographic': pv_projection = ccrs.Orthographic()
		if proj_name == 'Robinson':     pv_projection = ccrs.Robinson()
		if proj_name == 'Satellite':    pv_projection = ccrs.NearsidePerspective(central_longitude=17.5,central_latitude=36.4,satellite_height=6779000)
		if proj_name == 'Eckert IV':    pv_projection = ccrs.EckertIV()
		if proj_name == 'Equal Earth':  pv_projection = ccrs.EqualEarth()
		if proj_name == 'EPSG 3857':    pv_projection = epsg(3857)

		xyz_new = projection.transform_points(pv_projection,xyz[:,0], xyz[:,1])

		# Set minimum and maximum
		z_min = np.nanmin(data); z_max = np.nanmax(data)
		cbar_min = mcbar_min if mcbar_min >= -1e10 else z_min
		cbar_max = mcbar_max if mcbar_max <= 1e10  else z_max

		# Set extend
		extend = 'neither'
		if (cbar_min > z_min): extend = 'min'
		if (cbar_max < z_max): extend = 'max'
		if (cbar_min > z_min and cbar_max < z_max): extend = 'both'

		# Create triangulation
		tri = matplotlib.tri.Triangulation(xyz_new[:,0], xyz_new[:,1], triangles=triang)

		# Plot 
		if not mcbar_log:
			cplot = ax.tricontourf(tri,data,
						cmap=cmap,
						levels=np.linspace(max(cbar_min,z_min),min(cbar_max,z_max),256),
						norm=matplotlib.colors.Normalize(cbar_min,cbar_max),
						extend=extend,
						transform=projection
					)
		else:
			cplot = ax.tricontourf(tri,data,
						cmap=cmap,
						levels=np.linspace(max(cbar_min,z_min),min(cbar_max,z_max),256),
						norm=matplotlib.colors.LogNorm(cbar_min,cbar_max),
						transform=projection
					)

		# Colorbar
		if mdraw_colorbar:
			cbar = figure.colorbar(cplot,orientation='horizontal')
			cbar.set_label(label=mcbar_label,size=mcbar_font,weight='bold' if mcbar_bold else None,style='italic' if mcbar_ital else None)
			cbar.ax.tick_params(labelsize=mcbar_font)
			if not mcbar_log:
				cbar.locator = matplotlib.ticker.LinearLocator(numticks=10)#FixedLocator(np.linspace(cbar_min,cbar_max,10))
			else:
				cbar.locator = matplotlib.ticker.LogLocator()
			cbar.formatter = matplotlib.ticker.FormatStrFormatter('%.2f')
			cbar.update_ticks()

		# Save figure
		if msavefigure and mfilename:
			# Fix ugly white lines in PDF
			for c in cplot.collections: c.set_edgecolor("face")
			# Save figure
			figure.savefig(mfilename,dpi=moutdpi,bbox_inches='tight')

		return python_view.figure_to_image(figure)