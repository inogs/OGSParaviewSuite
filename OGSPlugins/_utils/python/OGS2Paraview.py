#!/usr/bin/env pvpython
#
# Generation of the OGS master file (.ogs) and the
# required mesh files to open simulation in ParaView.
#
# Usage: 
#	OGS2ParaView [-h] -n NAME -p PATH [--ave_freq_1 AVE_FREQ_1]
#                    [--gen-mesh] [-m MAP] [-l LIB]
#
# Arguments:
#		-n NAME, --name NAME    	Name of the simulation
#		-p PATH, --path PATH    	Path to the simulation files
#		--ave_freq_1 AVE_FREQ_1 	Read AVE_FREQ_1 instead of AVE_FREQ_2
#		--gen-mesh              	Force generation of the mesh file
#		-m MAP, --map MAP       	Projection type (default: merc)
#		-l LIB, --lib LIB       	Path to the libOGSmesh.so library
#
# Arnau Miro, OGS (2018)

import os, sys, datetime, argparse
import numpy as np, glob as gl, ConfigParser as cp

import OGSmesh

class OGS2ParaView:
	'''
	Class that handles writing the OGS main file. It passes from
	NetCDF (meshmask.nc) to the OGS master file (.ogs) and the OGS
	mesh file (.ogsmsh) that also includes the masks.
	'''
	def __init__(self,name,path,mode=1,meshmask="meshmask.nc",config=os.path.join(os.path.dirname(sys.argv[0]),"default.ini")):
		'''
		Class constructor.
		'''
		# Filenames
		self._name      = name
		self._path      = path
		self._meshmask  = meshmask
		self.masterfile = os.path.join(path,"%s.ogs" % name)
		self.meshfile   = os.path.join(path,"%s.ogsmsh" % name)
		self.cfgfile    = config
		# Mode
		self._mode = mode
		# Dictionaries
		self.vardicts = [dict() for i in xrange(5)]
		self.time_list = None
		# Open file for reading
		self.fid = open(self.masterfile,"w")

	def __del__(self):
		'''
		Class destructor.
		'''
		# Force close of the file
		self.fid.close()

	def checkExist(self,path,searchstr):
		'''
		Check if a path with a search string exists. In this case
		return the file, otherwise return None
		'''
		file = gl.glob(os.path.join(path,searchstr))
		if len(file) == 0:
			return None
		return os.path.basename(file[0])

	def readConfig(self):
		'''
		Read and load the configuration file for the variables dictionary
		if the user has specified one, otherwise load the default.
		'''
		if not os.path.exists(self.cfgfile):
			raise ValueError("Configuration file in <%s> not found!\nAborting..." % self.cfgfile)
		# Open config parser
		cfgpar = cp.ConfigParser()
		cfgpar.read(self.cfgfile)

		ii = 0
		for sect in ["AVE_PHYS","AVE_FREQ_1","AVE_FREQ_2","FORCINGS","GENERALS"]:
			self.vardicts[ii]["format"] = cfgpar.get(sect,"file format")
			self.vardicts[ii]["refvar"] = cfgpar.get(sect,"reference variable")
			self.vardicts[ii]["forcen"] = cfgpar.get(sect,"force var name") in ["true","True","yes","1"]
			self.vardicts[ii]["varlis"] = [var.split(',') for var in cfgpar.get(sect,"variable list").split('\n')]
			ii += 1

	def genMesh(self,force_mesh=False):
		'''
		Generate the binary .ogsmsh file.
		'''
		if not os.path.exists(self.meshfile) or force_mesh:
			mesh = OGSmesh.OGSmesh(self._path,maskname=self._meshmask)
			mesh.createOGSMesh("%s.ogsmsh" % self._name,self._path)

	def readTimeInstants(self):
		'''
		Read the time instants according to the operating mode.
		'''
		# Set index according to mode
		ind = 0
		if self._mode == 2: ind = 1
		if self._mode == 3: ind = 3
		if self._mode == 4: ind = 4
		if self.vardicts[ind]["refvar"] == '':
			raise ValueError("No reference value has been inputted!\nCheck your configuration file...")
		# Build the search string
		searchstr = self.vardicts[ind]["format"] % ("*",self.vardicts[ind]["refvar"])
		# Build the time list
		self.time_list = sorted([os.path.basename(file).split('.')[1] for file in gl.iglob(os.path.join(self._path,searchstr))])

	def listVars(self,idx):
		'''
		Lists all the items in a path given a search string with
		a wildcard using glob. Looks for the variable names.
		'''
		# Build the search string
		searchstr = os.path.join(self._path,self.vardicts[idx]["format"] % (self.time_list[0],"*"))
		# Return the variables list
		return sorted([os.path.basename(file).split('.')[2] for file in gl.iglob(os.path.join(self._path,searchstr))])

	def writeHeader(self):
		'''
		Header lets us know for which simulation and at which date 
		the file has been created.
		'''
		self.fid.write("# File for %s created at %s\n" % (self._name,datetime.datetime.now().strftime("%Y-%m-%d %H:%M")) )
		self.fid.write("# Generated using %s\n\n" % os.path.basename(sys.argv[0]))

	def writeWorkdir(self):
		'''
		Working directory sets the absolute path to where the simulation is
		stored. All subsequent paths are given as relative to the working
		directory path.
		'''
		self.fid.write("WRKDIR\n")
		self.fid.write("%s\n\n" % (os.path.abspath(self._path)))

	def writeMesh(self):
		'''
		Mesh sets the path to the .ogsmsh file, a binary file containing all the 
		mesh information and the basins and coasts mask.
		'''
		self.fid.write("MESH\n")
		self.fid.write("%s.ogsmsh : %s\n\n" % (self._name,self._meshmask))

	def writeVars(self,header,idx):
		'''
		Physical variables are written here. They are stored in the *.phys.* file and
		listed in a dictionary (ave_phys) at the beginning of this file.
		'''
		self.fid.write("%s\n" % header)
		if self.time_list == None:
			self.fid.write("0\n\n")
			return
		# List the variables
		varname = self.vardicts[idx]["varlis"] if self.vardicts[idx]["forcen"] else self.listVars(idx)
		cdfname = varname
		# Fix variable name and netcdf name in case force var name is active
		if self.vardicts[idx]["forcen"]:
			varname = [v[0].strip() for v in varname]
			cdfname = [v[1].strip() for v in cdfname]

		# Check if all the variables exist for the specified timesteps
		bwrite_var = np.any([os.path.exists(os.path.join(self._path,self.vardicts[idx]["format"] % (t,v))) for t in self.time_list for v in varname]) \
			if not self.vardicts[idx]["forcen"] \
			else np.any([os.path.exists(os.path.join(self._path,self.vardicts[idx]["format"] % (t,self.vardicts[idx]["refvar"]))) for t in self.time_list])

		# Write the variables
		if (len(self.time_list) == 0 or not bwrite_var):
			# If there are no instants or the instant or we couldn't find all the time instants
			# do not write the variable
			self.fid.write("0\n")
		else:
			# Write the number of variables
			self.fid.write("%d\n" % len(varname))
			# Write the variables
			for c,v in zip(varname,cdfname):
				fname = self.vardicts[idx]["format"] % ("*",v) if not self.vardicts[idx]["forcen"] else \
						self.vardicts[idx]["format"] % ("*",self.vardicts[idx]["refvar"])
				self.fid.write("%s : %s : %s\n" % (v,c,fname))
		self.fid.write("\n")

	def writeTime(self):
		'''
		Time instants are scanned and written here in the file.
		'''
		self.fid.write("TIME\n")
		if self.time_list == None:
			self.fid.write("0\n")
			return
		self.fid.write("%d\n" % len(self.time_list))
		for t in self.time_list:
			self.fid.write("%s\n" % t)
		self.fid.write("\n")

	def writeOGSFile(self,force_mesh=False):
		'''
		Writes OGS file according to the mode the user has inputted.
		'''
		# First write the header and the working directory
		self.writeHeader()
		self.writeWorkdir()
		
		# Then generate and write the mesh
		self.genMesh(True if self._mode == 0 else force_mesh)
		self.writeMesh()

		# Mode = 0 just generate the mesh
		self.readConfig()
		if self._mode > 0:
			self.readTimeInstants()		

		# Write the variables
		self.writeVars("AVE_PHYS",0)
		self.writeVars("AVE_FREQ",1 if self._mode == 2 else 2)
		self.writeVars("FORCINGS",3)
		self.writeVars("GENERALS",4)	

		# Write the timesteps
		self.writeTime() 



if __name__ == '__main__':
	'''
	Main implementation, with ARGPARSE, to generate the 
	OGS main file for ParaView.
	'''
	# Arguments
	argpar = argparse.ArgumentParser(prog="OGS2ParaView",description="Generation of OGS master file.")
	argpar.add_argument('-n','--name',type=str,help='Name of the simulation',required=True,dest='name')
	argpar.add_argument('-p','--path',type=str,help='Path to the simulation files',required=True,dest='path')
	argpar.add_argument('-m','--mode',type=int,help='Mode of operation',dest='mode')
	argpar.add_argument('-c','--config',type=str,help='Configuration file for the simulation',dest='conf')
	argpar.add_argument('--gen-mesh',action='store_true',help='Force generation of the mesh file',dest='gen_mesh')
	argpar.add_argument('--meshmask',type=str,help='Name of the meshmask to use',dest='mesh')

	# parse input arguments
	args=argpar.parse_args()
	if args.mode == None: args.mode = 1
	if not args.mesh: args.mesh = "meshmask.nc"
	if not args.conf: args.conf = os.path.join(os.path.dirname(sys.argv[0]),"default.ini")

	# Create an instance of the class
	OGS2P = OGS2ParaView(name=args.name,path=args.path,mode=args.mode,meshmask=args.mesh,config=args.conf)
	OGS2P.genMesh(force_mesh=args.gen_mesh)

	OGS2P.writeOGSFile()