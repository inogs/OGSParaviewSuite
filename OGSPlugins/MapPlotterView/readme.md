## Map Plotter View

The **Map Plotter View** is a python view designed to interface with bit.sea [Map Plotter](https://github.com/inogs/bit.sea/tree/master/MapPlotter) tool. This view can provide high resolution beautiful 2D maps directly from NetCDF files (without having to load into ParaView) using [cartopy](https://scitools.org.uk/cartopy/docs/latest/). The user can set all the parameters of the Map Plotter class through the ParaView interface.

This tool depends on:
* bit.sea [Map Plotter](https://github.com/inogs/bit.sea/tree/master/MapPlotter)
* the [PROJ library](https://proj.org/)
* [Cartopy](https://scitools.org.uk/cartopy/docs/latest/)
* the [requests](https://pypi.org/project/requests/) module
