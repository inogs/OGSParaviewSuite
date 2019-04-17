 pv_plugin(OGSSelectTools
 
   # Provide brief description for the plugin used as documentation for
   # the PARAVIEW_BUILD_PLUGIN_<PluginName> cmake option provided to the user.
   DESCRIPTION "Lets the user select basins from the basins mask"
 
   # If you want the plugin to be auto-loaded when ParaView starts, specify this option.
   # Users can manually mark any plugin to be auto-loaded using the Plugin Manager dialog.
   # This option is ignore for static-builds. All enabled plugins are auto-loaded in static
   # builds.
   AUTOLOAD
   
   # Specify this option if PARAVIEW_BUILD_PLUGIN_<PluginName> option should default to TRUE.
   # If not specified, it defaults to FALSE and the user must turn it ON to build this plugin.
   # Note the user can always turn PARAVIEW_BUILD_PLUGIN_<PluginName> off using cmake.
   DEFAULT_ENABLED

   # If providing more than 1 plugin or plugin is named differently (in add_paraview_plugin call)
   # than the <PluginName> specified,
   # you can use this option to notify ParaView of the plugin library names. ParaView uses these
   # names to locate the plugin at run time.
   #PLUGIN_NAMES OGSSelectCoast OGSSelectBasin
 )