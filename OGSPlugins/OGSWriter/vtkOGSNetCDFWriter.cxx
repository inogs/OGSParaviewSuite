/*=========================================================================

	Program:   OGSWriter
	Module:    vtkOGSNetCDFWriter.cxx

	Copyright (c) 2019 Arnau Miro, OGS
	All rights reserved.

		 This software is distributed WITHOUT ANY WARRANTY; without even
		 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
		 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkOGSNetCDFWriter.h"

#include "vtkCell.h"
#include "vtkPoints.h"
#include "vtkCallbackCommand.h"
#include "vtkCommand.h"
#include "vtkExecutive.h"
#include "vtkCellData.h"
#include "vtkPointData.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkTypeUInt8Array.h"
#include "vtkFloatArray.h"
#include "vtkAbstractArray.h"
#include "vtkDoubleArray.h"
#include "vtkStringArray.h"
#include "vtkCellLocator.h"
#include "vtkGenericCell.h"
#include "vtkRectilinearGrid.h"
#include "vtkUnstructuredGrid.h"
#include "vtkPolyData.h"

#include "vtkObjectFactory.h"

#include<vector>
#include<string>

#define CELL_TOLERANCE_FACTOR_SQR 1e-6

vtkStandardNewMacro(vtkOGSNetCDFWriter);

//----------------------------------------------------------------------------

/*
	Macro to set the array precision 
*/
#define FLDARRAY double
#define FLDMASK  uint8_t
#define VTKARRAY vtkDoubleArray
#define VTKMASK  vtkTypeUInt8Array

#include "../_utils/V3.h"
#include "../_utils/field.h"
#include "../_utils/vtkFields.hpp"
#include "../_utils/vtkOperations.hpp"
#include "../_utils/OGS.hpp"
#include "../_utils/netcdfio.hpp"

//----------------------------------------------------------------------------
void FindCell(v3::V3 p, v3::V3v &xyz, int &pos_min, double &dist_min, double epsi) {

	int ii_st = pos_min + 1, ii_ed = pos_min+xyz.len()/100;
	if (ii_ed > xyz.len()) ii_ed = xyz.len();

	// Do a reduced serch on xyz and compute the distance to p
	for (int ii = ii_st; ii < ii_ed; ++ii) {
		v3::V3 dist = p - xyz[ii];
		double d = dist.norm2();
		if (d < dist_min) { dist_min = d; pos_min = ii; }
		if (dist_min < epsi) return; // Exit if we are under tolerance
	}
	// If we don't fit the tolerance expectations, loop the rest of the mesh
	// and try to find the value
	for(int ii = ii_ed; ii < xyz.len(); ++ii) {
		v3::V3 dist = p - xyz[ii];
		double d = dist.norm2();
		if (d < dist_min) { dist_min = d; pos_min = ii; }
		if (dist_min < epsi) return; // Exit if we are under tolerance
	}
	// If we don't fit the tolerance expectations, loop the rest of the mesh
	// and try to find the value
	for(int ii = 0; ii < ii_st; ++ii) {
		v3::V3 dist = p - xyz[ii];
		double d = dist.norm2();
		if (d < dist_min) { dist_min = d; pos_min = ii; }
		if (dist_min < epsi) return; // Exit if we are under tolerance
	}
	// Assume that we have found the point at this point
}

void FindPoint(v3::V3 p, v3::V3v &xyz, int &pos_min, double &dist_min, double epsi) {

	int ii_st = pos_min - xyz.len()/100, ii_ed = pos_min+xyz.len()/100;
	if (ii_st < 0)         ii_st = 0;
	if (ii_ed > xyz.len()) ii_ed = xyz.len();

	// Do a reduced serch on xyz and compute the distance to p
	for (int ii = ii_st; ii < ii_ed; ++ii) {
		v3::V3 dist = p - xyz[ii];
		double d = dist.norm2();
		if (d < dist_min) { dist_min = d; pos_min = ii; }
		if (dist_min < epsi) return; // Exit if we are under tolerance
	}
	// If we don't fit the tolerance expectations, loop the rest of the mesh
	// and try to find the value
	for(int ii = ii_ed; ii < xyz.len(); ++ii) {
		v3::V3 dist = p - xyz[ii];
		double d = dist.norm2();
		if (d < dist_min) { dist_min = d; pos_min = ii; }
		if (dist_min < epsi) return; // Exit if we are under tolerance
	}
	// If we don't fit the tolerance expectations, loop the rest of the mesh
	// and try to find the value
	for(int ii = 0; ii < ii_st; ++ii) {
		v3::V3 dist = p - xyz[ii];
		double d = dist.norm2();
		if (d < dist_min) { dist_min = d; pos_min = ii; }
		if (dist_min < epsi) return; // Exit if we are under tolerance
	}
	// Assume that we have found the point at this point
}

//----------------------------------------------------------------------------
int InterpolateData(vtkDataSet *input, vtkDataSet *output, double dfact, bool iscelld) {
	// Interpolate from input to output for each point of the input, 
	// find its equivalent on the output

	// Loop the number of points in the input mesh
	double epsi = 1.e-6;
	vtkAbstractArray *inArray, *outArray;

	if (iscelld) {
		// Compute the cell centers
		v3::V3v xyz1 = VTK::getVTKCellCenters(input,dfact);
		v3::V3v xyz2 = VTK::getVTKCellCenters(output,dfact);
		
		int cellId = -1;
		int nvars = output->GetCellData()->GetNumberOfArrays();

		for (int ii = 0; ii < xyz1.len(); ++ii) {
			double dist_min = 999.;

			FindCell(xyz1[ii],xyz2,cellId,dist_min,epsi);
			if (epsi < dist_min) epsi = dist_min;

			// Interpolate data from input to output
			for (int varId = 0; varId < nvars; ++varId) {
				// Recover output array
				outArray = output->GetCellData()->GetArray(varId);
				std::string outName = std::string(outArray->GetName());
				// Recover input array
				inArray = input->GetCellData()->GetArray(outName.c_str());
				// Interpolate data
				output->GetCellData()->CopyTuple(inArray,outArray,ii,cellId);
			}
		}
	} else {
		// Compute the cell centers
		v3::V3v xyz1 = VTK::getVTKCellPoints(input,dfact);
		v3::V3v xyz2 = VTK::getVTKCellPoints(output,dfact);
	
		int cId = 0, cellId[2] = {-1,-1};
		int nvars = output->GetPointData()->GetNumberOfArrays();

		for (int ii = 0; ii < xyz1.len(); ++ii) {
			double dist_min = 999.;

			if (ii > 0 && std::fabs(xyz1[ii-1][2]-xyz1[ii][2]) > 1e-3)
				cId = (cId == 0) ? 1 : 0;

			FindPoint(xyz1[ii],xyz2,cellId[cId],dist_min,epsi);
			if (epsi < dist_min) epsi = dist_min;

			// Interpolate data from input to output
			for (int varId = 0; varId < nvars; ++varId) {
				// Recover output array
				outArray = output->GetPointData()->GetArray(varId);
				std::string outName = std::string(outArray->GetName());
				// Recover input array
				inArray = input->GetPointData()->GetArray(outName.c_str());
				// Interpolate data
				output->GetPointData()->CopyTuple(inArray,outArray,ii,cellId[cId]);
			}		
		}
	}
	return 1;
}

//----------------------------------------------------------------------------
vtkOGSNetCDFWriter::vtkOGSNetCDFWriter() : FileName(nullptr), OGSFile(nullptr), varname(nullptr), dfact(1000.), singlevar(0), SaveAll(0) {}

//----------------------------------------------------------------------------
vtkOGSNetCDFWriter::~vtkOGSNetCDFWriter() {
	this->SetFileName(nullptr);
	this->SetOGSFile(nullptr);
	this->Setvarname(nullptr);
}

//----------------------------------------------------------------------------
int vtkOGSNetCDFWriter::FillInputPortInformation( int vtkNotUsed(port), vtkInformation* info) {
	info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
	return 1;
}

//----------------------------------------------------------------------------
void vtkOGSNetCDFWriter::WriteData() {
	// Make sure we only export from one mesh
	if(this->GetNumberOfInputConnections(0) != 1)
		vtkErrorMacro("Exactly one input required.");

	// Recover the input
	vtkDataSet *input = vtkDataSet::SafeDownCast( this->GetExecutive()->GetInputData(0,0) );

	// Check that varname is not metadata
	if (std::string(this->varname) == std::string("Metadata")) {
		vtkWarningMacro("Metadata cannot be the input variable! Aborting.");
		return;
	}

	// Recover Metadata array (depth factor)
	vtkStringArray *vtkmetadata = vtkStringArray::SafeDownCast(
		input->GetFieldData()->GetAbstractArray("Metadata"));
	this->dfact = (vtkmetadata != NULL) ? std::stod( vtkmetadata->GetValue(2) ) : this->dfact;
	std::string aux = (vtkmetadata != NULL) ? vtkmetadata->GetValue(4) : "";
	if (aux != std::string("")) { this->OGSFile = new char[aux.length()+1]; std::strcpy(this->OGSFile,aux.c_str()); }
	this->projId = (vtkmetadata != NULL) ? std::stod( vtkmetadata->GetValue(7) ) : 0;
		
	if (vtkmetadata == NULL) 
		vtkWarningMacro("Field array Metadata not found! Using user input value.");

	// Try to load as Cell data
	// Use the variable in varname to check if we are using CellData or PointData
	bool iscelld = true; vtkAbstractArray *vtkArray;
	vtkArray = input->GetCellData()->GetArray(this->varname);
	// Then array might be point data
	if (!vtkArray) {
		vtkArray = input->GetPointData()->GetArray(this->varname);
		iscelld = false;
	}

	// Deal with the kind of data
	if ( input->IsA("vtkRectilinearGrid") ) {
		// Rectilinear Grid
		if (this->singlevar)
			this->writeNetCDFRectilinearGrid1Var(vtkRectilinearGrid::SafeDownCast(input),vtkArray,iscelld);
		else
			this->writeNetCDFRectilinearGridnVar(vtkRectilinearGrid::SafeDownCast(input),iscelld);
	} else if ( input->IsA("vtkUnstructuredGrid") ) {
		// Unstructured Grid
		if (this->singlevar)
			this->writeNetCDFUnstructuredGrid1Var(vtkUnstructuredGrid::SafeDownCast(input),vtkArray,iscelld);
		else
			this->writeNetCDFUnstructuredGridnVar(vtkUnstructuredGrid::SafeDownCast(input),iscelld);
	} else if ( input->IsA("vtkPolyData") ) {
		// Poly Data
		if (this->singlevar)
			this->writeNetCDFPolyData1Var(vtkPolyData::SafeDownCast(input),vtkArray,iscelld);
		else
			this->writeNetCDFPolyDatanVar(vtkPolyData::SafeDownCast(input),iscelld);
	} else {
		vtkWarningMacro("Type of data not recognized to generate a NetCDF file!");
	}
}

//----------------------------------------------------------------------------
void vtkOGSNetCDFWriter::writeNetCDFRectilinearGrid1Var(vtkRectilinearGrid *input, vtkAbstractArray *vtkArray, bool iscelld) {
	// Recover mesh dimensions
	int nx = input->GetXCoordinates()->GetNumberOfTuples();
	int ny = input->GetYCoordinates()->GetNumberOfTuples();
	int nz = input->GetZCoordinates()->GetNumberOfTuples();
	
	// Compute the cell centers
	v3::V3v xyz = (iscelld) ? VTK::getVTKCellCenters(input,this->dfact) : VTK::getVTKCellPoints(input,this->dfact);

	// Load the data according to the variable name. Data will be converted to float
	// and stored to the NetCDF file
	field::Field<float> data;

	// Load the variable
	if (std::string(this->varname) == std::string("basins mask")      || 
		std::string(this->varname) == std::string("coast mask")       ||
		std::string(this->varname) == std::string("Okubo-Weiss mask") ||
		std::string(this->varname) == std::string("Q-criterion mask")) {
		// Use mask uint8
		field::Field<FLDMASK> array;
		array = VTK::createFieldfromVTK<VTKMASK,FLDMASK>( VTKMASK::SafeDownCast(vtkArray) );
		data = array.convert<float>();
	} else {
		// Use normal array
		field::Field<FLDARRAY> array;
		array = VTK::createFieldfromVTK<VTKARRAY,FLDARRAY>( VTKARRAY::SafeDownCast(vtkArray) );
		data = array.convert<float>();
	}

	// Write NetCDF file
	int retval;
	if (iscelld) {
		int dims[] = {nx-1,ny-1,nz-1};
		retval = NetCDF::writeNetCDF(this->FileName,this->varname,dims,xyz,data);
	} else {
		int dims[] = {nx,ny,nz};
		retval = NetCDF::writeNetCDF(this->FileName,this->varname,dims,xyz,data);
	}
	if ( retval != NETCDF_OK ) vtkWarningMacro("Unexpected error writing NetCDF file!");
}

//----------------------------------------------------------------------------
void vtkOGSNetCDFWriter::writeNetCDFRectilinearGridnVar(vtkRectilinearGrid *input, bool iscelld) {
	// Recover mesh dimensions
	int nx = input->GetXCoordinates()->GetNumberOfTuples();
	int ny = input->GetYCoordinates()->GetNumberOfTuples();
	int nz = input->GetZCoordinates()->GetNumberOfTuples();
	
	// Compute the cell centers
	v3::V3v xyz = (iscelld) ? VTK::getVTKCellCenters(input,this->dfact) : VTK::getVTKCellPoints(input,this->dfact);

	// Do a first scan of the loaded variables and create a vector 
	// containing the variable names
	vtkAbstractArray *vtkArray;
	std::vector<std::string> varnames;
	std::vector<int> varId;

	int nvars = (iscelld) ? input->GetCellData()->GetNumberOfArrays() : input->GetPointData()->GetNumberOfArrays();
	
	for (int vid = 0; vid < nvars; ++vid) {
		// Retrieve the variable name
		vtkArray = (iscelld) ? input->GetCellData()->GetArray(vid) : input->GetPointData()->GetArray(vid);
		std::string vname = std::string(vtkArray->GetName());
		// Discard some variables
		if ( !this->SaveAll && vname == std::string("basins mask")      ) continue;
		if ( !this->SaveAll && vname == std::string("coast mask")       ) continue;
		if ( !this->SaveAll && vname == std::string("Okubo-Weiss mask") ) continue;
		if ( !this->SaveAll && vname == std::string("Q-criterion mask") ) continue;
		if ( !this->SaveAll && vname == std::string("coast mask")       ) continue;
		if ( !this->SaveAll && vname == std::string("e1")               ) continue;
		if ( !this->SaveAll && vname == std::string("e2")               ) continue;
		if ( !this->SaveAll && vname == std::string("e3")               ) continue;
		// Store the variable name in the vector array
		int ncomp = vtkArray->GetNumberOfValues()/vtkArray->GetNumberOfTuples();
		if (ncomp > 1) {
			for (int cid = 0; cid < ncomp; ++cid) {
				std::string auxname = vname + std::string("_") + std::to_string(cid);
				std::replace(auxname.begin(),auxname.end(),' ','_'); // Replace white spaces
				varnames.push_back( auxname );
				varId.push_back( vid );
			}
		} else {
			varnames.push_back( vname );
			varId.push_back( vid );
		}
	}

	// Now that we have the variables to store, let's pack them into a single
	// vector of float arrays to later save the NetCDF file
	nvars = (int)(varnames.size());
	std::vector<float*> data(nvars);

	for (int vid = 0; vid < nvars; ++vid) {
		// Retrieve the variable
		vtkArray = (iscelld) ? input->GetCellData()->GetArray(varId[vid]) : input->GetPointData()->GetArray(varId[vid]);
		std::string vname = std::string(vtkArray->GetName());
		field::Field<float> aux;
		// Load the variable
		if (vname == std::string("basins mask") || vname == std::string("coast mask") ||
			vname == std::string("Okubo-Weiss mask") || vname == std::string("Q-criterion mask")) {
			// Use mask uint8
			field::Field<FLDMASK> array;
			array = VTK::createFieldfromVTK<VTKMASK,FLDMASK>( VTKMASK::SafeDownCast(vtkArray) );
			aux = array.convert<float>();
		} else {
			// Use normal array
			field::Field<FLDARRAY> array;
			array = VTK::createFieldfromVTK<VTKARRAY,FLDARRAY>( VTKARRAY::SafeDownCast(vtkArray) );
			aux = array.convert<float>();
		}
		// Copy the array to data
		for (int jj = 0; jj < aux.get_m(); ++jj) {
			data[vid + jj] = new float[aux.get_n()];
			for (int ii = 0; ii < aux.get_n(); ++ii)
				data[vid + jj][ii] = aux[ii][jj];
		}
		// Actualize index in case of multidimensional array
		vid += aux.get_m() - 1;
	}

	// Write NetCDF file
	int retval;
	if (iscelld) {
		int dims[] = {nx-1,ny-1,nz-1};
		retval = NetCDF::writeNetCDF(this->FileName,varnames.data(),nvars,dims,xyz,data.data());
	} else {
		int dims[] = {nx,ny,nz};
		retval = NetCDF::writeNetCDF(this->FileName,varnames.data(),nvars,dims,xyz,data.data());
	}
	if ( retval != NETCDF_OK ) vtkWarningMacro("Unexpected error writing NetCDF file!");
}

//----------------------------------------------------------------------------
void vtkOGSNetCDFWriter::writeNetCDFUnstructuredGrid1Var(vtkUnstructuredGrid *input, vtkAbstractArray *vtkArray, bool iscelld) {
	// Since the NetCDF files can only store structured data, the rectilinear mesh should
	// be recreated at this point. We can do so using the metadata array and the OGS class.
	ogs::OGS ogsdata;
	ogsdata.SetFile(this->OGSFile);
	if (ogsdata.readMainFile() == -1) { 
		vtkWarningMacro("Cannot read <"<<this->OGSFile<<">! Aborting.");
		return;
	}

	// Now read the mesh file
	if (ogsdata.readMesh(this->projId) < 0) {
		vtkWarningMacro("Problems reading the mesh file! Aborting.");
		return;
	}

	// Create the auxiliar rectilinear grid
	vtkRectilinearGrid *rgrid_mesh;
	rgrid_mesh = vtkRectilinearGrid::New();

	VTK::createRectilinearGrid(ogsdata.nlon(),
							   ogsdata.nlat(),
							   ogsdata.nlev(),
							   ogsdata.lon2meters(),
							   ogsdata.lat2meters(),
							   ogsdata.nav_lev(),
							   this->dfact,
							   rgrid_mesh
							  );

	// Generate an array for the rectilinear grid
	int n = (iscelld) ? ogsdata.ncells() : ogsdata.nlon()*ogsdata.nlat()*ogsdata.nlev();
	int ncomp = vtkArray->GetNumberOfValues()/vtkArray->GetNumberOfTuples();
	std::string vname = std::string( vtkArray->GetName() );

	// Load the variable
	if (vname == std::string("basins mask") || vname == std::string("coast mask") ||
		vname == std::string("Okubo-Weiss mask") || vname == std::string("Q-criterion mask")) {
		// Use mask uint8
		FLDMASK aux = 0;
		field::Field<FLDMASK> array(n,ncomp,aux);
		VTKMASK *vtkFArray; vtkFArray = VTK::createVTKfromField<VTKMASK,FLDMASK>(vname.c_str(),array);
		if (iscelld)
			rgrid_mesh->GetCellData()->AddArray(vtkFArray);
		else
			rgrid_mesh->GetPointData()->AddArray(vtkFArray);
	} else {
		// Use normal array
		field::Field<FLDARRAY> array(n,ncomp,MISSING_VALUE);
		VTKARRAY *vtkFArray; vtkFArray = VTK::createVTKfromField<VTKARRAY,FLDARRAY>(vname.c_str(),array);
		if (iscelld)
			rgrid_mesh->GetCellData()->AddArray(vtkFArray);
		else
			rgrid_mesh->GetPointData()->AddArray(vtkFArray);
	}

	// Now interpolate the data in the new rectilinear grid
	InterpolateData(input,rgrid_mesh,this->dfact,iscelld);

	// Use the rectilinear grid save algorithm to save the NetCDF file
	vtkArray = (iscelld) ? rgrid_mesh->GetCellData()->GetArray(vname.c_str()) : rgrid_mesh->GetPointData()->GetArray(vname.c_str());
	this->writeNetCDFRectilinearGrid1Var(rgrid_mesh,vtkArray,iscelld);

	// Free some memory
	rgrid_mesh->Delete();
}

//----------------------------------------------------------------------------
void vtkOGSNetCDFWriter::writeNetCDFUnstructuredGridnVar(vtkUnstructuredGrid *input, bool iscelld) {
	// Since the NetCDF files can only store structured data, the rectilinear mesh should
	// be recreated at this point. We can do so using the metadata array and the OGS class.
	ogs::OGS ogsdata;
	ogsdata.SetFile(this->OGSFile);
	if (ogsdata.readMainFile() == -1) { 
		vtkWarningMacro("Cannot read <"<<this->OGSFile<<">! Aborting.");
		return;
	}

	// Now read the mesh file
	if (ogsdata.readMesh(this->projId) < 0) {
		vtkWarningMacro("Problems reading the mesh file! Aborting.");
		return;
	}

	// Create the auxiliar rectilinear grid
	vtkRectilinearGrid *rgrid_mesh;
	rgrid_mesh = vtkRectilinearGrid::New();

	VTK::createRectilinearGrid(ogsdata.nlon(),
							   ogsdata.nlat(),
							   ogsdata.nlev(),
							   ogsdata.lon2meters(),
							   ogsdata.lat2meters(),
							   ogsdata.nav_lev(),
							   this->dfact,
							   rgrid_mesh
							  );

	// Generate the arrays for the rectilinear grid
	int n = (iscelld) ? ogsdata.ncells() : ogsdata.nlon()*ogsdata.nlat()*ogsdata.nlev();
	int nvars = (iscelld) ? input->GetCellData()->GetNumberOfArrays() : input->GetPointData()->GetNumberOfArrays();
	vtkAbstractArray *vtkArray;
	
	for (int vid = 0; vid < nvars; ++vid) {
		// Retrieve the variable name
		vtkArray = (iscelld) ? input->GetCellData()->GetArray(vid) : input->GetPointData()->GetArray(vid);
		std::string vname = std::string(vtkArray->GetName());
		int ncomp = vtkArray->GetNumberOfValues()/vtkArray->GetNumberOfTuples();

		// Load the variable
		if (vname == std::string("basins mask") || vname == std::string("coast mask") ||
			vname == std::string("Okubo-Weiss mask") || vname == std::string("Q-criterion mask")) {
			// Use mask uint8
			FLDMASK aux = 0;
			field::Field<FLDMASK> array(n,ncomp,aux);
			VTKMASK *vtkFArray; vtkFArray = VTK::createVTKfromField<VTKMASK,FLDMASK>(vname.c_str(),array);
			if (iscelld)
				rgrid_mesh->GetCellData()->AddArray(vtkFArray);
			else
				rgrid_mesh->GetPointData()->AddArray(vtkFArray);
		} else {
			// Use normal array
			field::Field<FLDARRAY> array(n,ncomp,MISSING_VALUE);
			VTKARRAY *vtkFArray; vtkFArray = VTK::createVTKfromField<VTKARRAY,FLDARRAY>(vname.c_str(),array);
			if (iscelld)
				rgrid_mesh->GetCellData()->AddArray(vtkFArray);
			else
				rgrid_mesh->GetPointData()->AddArray(vtkFArray);
		}
	}

	// Now interpolate the data in the new rectilinear grid
	InterpolateData(input,rgrid_mesh,this->dfact,iscelld);

	// Use the rectilinear grid save algorithm to save the NetCDF file
	this->writeNetCDFRectilinearGridnVar(rgrid_mesh,iscelld);

	// Free some memory
	rgrid_mesh->Delete();
}

//----------------------------------------------------------------------------
void vtkOGSNetCDFWriter::writeNetCDFPolyData1Var(vtkPolyData *input, vtkAbstractArray *vtkArray, bool iscelld) {
	if (iscelld) {
		vtkWarningMacro("vtkPolyData with cell data! Aborting...");
		return;
	}

	// Compute points
	v3::V3v xyz = VTK::getVTKCellPoints(input,this->dfact);

	// Load the data according to the variable name. Data will be converted to float
	// and stored to the NetCDF file
	field::Field<float> data;

	// Load the variable
	if (std::string(this->varname) == std::string("basins mask")      || 
		std::string(this->varname) == std::string("coast mask")       ||
		std::string(this->varname) == std::string("Okubo-Weiss mask") ||
		std::string(this->varname) == std::string("Q-criterion mask")) {
		// Use mask uint8
		field::Field<FLDMASK> array;
		array = VTK::createFieldfromVTK<VTKMASK,FLDMASK>( VTKMASK::SafeDownCast(vtkArray) );
		data = array.convert<float>();
	} else {
		// Use normal array
		field::Field<FLDARRAY> array;
		array = VTK::createFieldfromVTK<VTKARRAY,FLDARRAY>( VTKARRAY::SafeDownCast(vtkArray) );
		data = array.convert<float>();
	}

	// Write NetCDF file
	if ( NetCDF::writeNetCDFProfile(this->FileName,this->varname,xyz.len(),xyz,data) != NETCDF_OK ) 
		vtkWarningMacro("Unexpected error writing NetCDF file!");
}

//----------------------------------------------------------------------------
void vtkOGSNetCDFWriter::writeNetCDFPolyDatanVar(vtkPolyData *input, bool iscelld) {
	if (iscelld) {
		vtkWarningMacro("vtkPolyData with cell data! Aborting...");
		return;
	}

	// Compute points
	v3::V3v xyz = VTK::getVTKCellPoints(input,this->dfact);

	// Do a first scan of the loaded variables and create a vector 
	// containing the variable names
	vtkAbstractArray *vtkArray;
	std::vector<std::string> varnames;
	std::vector<int> varId;

	int nvars = (iscelld) ? input->GetCellData()->GetNumberOfArrays() : input->GetPointData()->GetNumberOfArrays();
	
	for (int vid = 0; vid < nvars; ++vid) {
		// Retrieve the variable name
		vtkArray = (iscelld) ? input->GetCellData()->GetArray(vid) : input->GetPointData()->GetArray(vid);
		std::string vname = std::string(vtkArray->GetName());
		// Discard some variables
		if ( !this->SaveAll && vname == std::string("basins mask")      ) continue;
		if ( !this->SaveAll && vname == std::string("coast mask")       ) continue;
		if ( !this->SaveAll && vname == std::string("Okubo-Weiss mask") ) continue;
		if ( !this->SaveAll && vname == std::string("Q-criterion mask") ) continue;
		if ( !this->SaveAll && vname == std::string("coast mask")       ) continue;
		if ( !this->SaveAll && vname == std::string("e1")               ) continue;
		if ( !this->SaveAll && vname == std::string("e2")               ) continue;
		if ( !this->SaveAll && vname == std::string("e3")               ) continue;
		// Store the variable name in the vector array
		int ncomp = vtkArray->GetNumberOfValues()/vtkArray->GetNumberOfTuples();
		if (ncomp > 1) {
			for (int cid = 0; cid < ncomp; ++cid) {
				std::string auxname = vname + std::string("_") + std::to_string(cid);
				std::replace(auxname.begin(),auxname.end(),' ','_'); // Replace white spaces
				varnames.push_back( auxname );
				varId.push_back( vid );
			}
		} else {
			varnames.push_back( vname );
			varId.push_back( vid );
		}
	}

	// Now that we have the variables to store, let's pack them into a single
	// vector of float arrays to later save the NetCDF file
	nvars = (int)(varnames.size());
	std::vector<float*> data(nvars);

	for (int vid = 0; vid < nvars; ++vid) {
		// Retrieve the variable
		vtkArray = (iscelld) ? input->GetCellData()->GetArray(varId[vid]) : input->GetPointData()->GetArray(varId[vid]);
		std::string vname = std::string(vtkArray->GetName());
		field::Field<float> aux;
		// Load the variable
		if (vname == std::string("basins mask") || vname == std::string("coast mask") ||
			vname == std::string("Okubo-Weiss mask") || vname == std::string("Q-criterion mask")) {
			// Use mask uint8
			field::Field<FLDMASK> array;
			array = VTK::createFieldfromVTK<VTKMASK,FLDMASK>( VTKMASK::SafeDownCast(vtkArray) );
			aux = array.convert<float>();
		} else {
			// Use normal array
			field::Field<FLDARRAY> array;
			array = VTK::createFieldfromVTK<VTKARRAY,FLDARRAY>( VTKARRAY::SafeDownCast(vtkArray) );
			aux = array.convert<float>();
		}
		// Copy the array to data
		for (int jj = 0; jj < aux.get_m(); ++jj) {
			data[vid + jj] = new float[aux.get_n()];
			for (int ii = 0; ii < aux.get_n(); ++ii)
				data[vid + jj][ii] = aux[ii][jj];
		}
		// Actualize index in case of multidimensional array
		vid += aux.get_m() - 1;
	}

	// Write NetCDF file
	if ( NetCDF::writeNetCDFProfile(this->FileName,varnames.data(),nvars,xyz.len(),xyz,data.data()) != NETCDF_OK ) 
		vtkWarningMacro("Unexpected error writing NetCDF file!");
}

//----------------------------------------------------------------------------
int vtkOGSNetCDFWriter::Write() {
	return Superclass::Write();
}