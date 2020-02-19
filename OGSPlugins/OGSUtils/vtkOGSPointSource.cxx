/*=========================================================================

	Program:   OGSUtils
	Module:    vtkOGSPointSource.cxx

	Copyright (c) 2018 Arnau Miro, OGS
	All rights reserved.

		 This software is distributed WITHOUT ANY WARRANTY; without even
		 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
		 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkOGSPointSource.h"

#include "vtkCellArray.h"
#include "vtkMath.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkRandomSequence.h"

#include "projection.h"

#include <cfloat>
#include <cmath>

vtkStandardNewMacro(vtkOGSPointSource);

//---------------------------------------------------------------------------
// Specify a random sequence, or use the non-threadsafe one in vtkMath by
// default.
vtkCxxSetObjectMacro(vtkOGSPointSource, RandomSequence, vtkRandomSequence);

//----------------------------------------------------------------------------
vtkOGSPointSource::vtkOGSPointSource(vtkIdType numPts) {
	this->NumberOfPoints = (numPts > 0 ? numPts : 1);

	this->Center[0] = 0.0;
	this->Center[1] = 0.0;
	this->Center[2] = 0.0;

	this->Radius = 0.;
	this->DepthScale = 1000.;
	this->Projection = 0;

	this->Distribution = VTK_POINT_UNIFORM;
	this->OutputPointsPrecision = SINGLE_PRECISION;
	this->RandomSequence = nullptr;

	this->SetNumberOfInputPorts(0);
}

//----------------------------------------------------------------------------
vtkOGSPointSource::~vtkOGSPointSource() {
	this->SetRandomSequence(nullptr);
}

//----------------------------------------------------------------------------
int vtkOGSPointSource::RequestData(vtkInformation *vtkNotUsed(request),
	vtkInformationVector **vtkNotUsed(inputVector), vtkInformationVector *outputVector) {

	// get the info object
	vtkInformation *outInfo = outputVector->GetInformationObject(0);

	// get the output
	vtkPolyData *output = vtkPolyData::SafeDownCast(
		outInfo->Get(vtkDataObject::DATA_OBJECT()));

	vtkIdType i;
	double theta, rho, cosphi, sinphi, radius;
	double x[3];
	vtkPoints *newPoints;
	vtkCellArray *newVerts;

	newPoints = vtkPoints::New();

	// Set the desired precision for the points in the output.
	if(this->OutputPointsPrecision == vtkAlgorithm::DOUBLE_PRECISION)
		newPoints->SetDataType(VTK_DOUBLE);
	else
		newPoints->SetDataType(VTK_FLOAT);

	newPoints->Allocate(this->NumberOfPoints);
	newVerts = vtkCellArray::New();
	newVerts->Allocate(newVerts->EstimateSize(1,this->NumberOfPoints));

	newVerts->InsertNextCell(this->NumberOfPoints);

	if (this->Distribution == VTK_POINT_SHELL) {  // only produce points on the surface of the sphere
		for (i=0; i<this->NumberOfPoints; i++) {
			cosphi = 1 - 2 * this->Random();
			sinphi = sqrt(1 - cosphi*cosphi);
			radius = this->Radius * sinphi;
			theta = 2.0 * vtkMath::Pi() * this->Random();
			x[0] = this->Center[0] + radius*cos(theta);
			x[1] = this->Center[1] + radius*sin(theta);
			x[2] = this->Center[2] + this->Radius*cosphi;
			newVerts->InsertCellPoint(newPoints->InsertNextPoint(x));
		}
	} else { // uniform distribution throughout the sphere volume
		for (i=0; i<this->NumberOfPoints; i++) {
			cosphi = 1 - 2*this->Random();
			sinphi = sqrt(1 - cosphi*cosphi);
			rho = this->Radius*pow(this->Random(),0.33333333);
			radius = rho * sinphi;
			theta = 2.0 * vtkMath::Pi() * this->Random();
			x[0] = this->Center[0] + radius*cos(theta);
			x[1] = this->Center[1] + radius*sin(theta);
			x[2] = this->Center[2] + rho*cosphi;
			newVerts->InsertCellPoint(newPoints->InsertNextPoint(x));
		}
	}

	// Update ourselves and release memory
	output->SetPoints(newPoints);
	newPoints->Delete();

	output->SetVerts(newVerts);
	newVerts->Delete();

	return 1;
}

//----------------------------------------------------------------------------
double vtkOGSPointSource::Random() {
	if (!this->RandomSequence)
		return vtkMath::Random();

	this->RandomSequence->Next();
	return this->RandomSequence->GetValue();
}

// ----------------------------------------------------------------------
void vtkOGSPointSource::SetLonLat(double lon, double lat) {
	// Conversion of lon, lat to a projected point
	PROJ::Projection p;
	switch(this->Projection) {
		case 0:
			p.transform_point("degrees","mercator",lon,lat);
			break;
		case 1:
			p.transform_point("degrees","cylindrical",lon,lat);
			break;
		case 2:
			p.transform_point("degrees","google",lon,lat);
			break;
		case 3:
			p.transform_point("degrees","mollweide",lon,lat);
			break;
		case 4:
			p.transform_point("degrees","orthographic",lon,lat);
			break;
		case 5:
			p.transform_point("degrees","robinson",lon,lat);
			break;
		case 6:
			p.transform_point("degrees","satellite",lon,lat);
			break;
		case 7:
			p.transform_point("degrees","eckert iv",lon,lat);
			break;
		case 8:
			p.transform_point("degrees","equal earth",lon,lat);
			break;
		case 9:
			p.transform_point("degrees","epsg 3857",lon,lat);
			break; 
	}
	// Set center
	this->Center[0] = lon;
	this->Center[1] = lat;
	this->Modified();
}
void vtkOGSPointSource::GetLonLat(double &lon, double &lat) {
	// Conversion to lon, lat
	PROJ::Projection p;
	lon = this->Center[0]; lat = this->Center[2];
	switch(this->Projection) {
		case 0:
			p.transform_point("mercator","degrees",lon,lat);
			break;
		case 1:
			p.transform_point("cylindrical","degrees",lon,lat);
			break;
		case 2:
			p.transform_point("google","degrees",lon,lat);
			break;
		case 3:
			p.transform_point("mollweide","degrees",lon,lat);
			break;
		case 4:
			p.transform_point("orthographic","degrees",lon,lat);
			break;
		case 5:
			p.transform_point("robinson","degrees",lon,lat);
			break;
		case 6:
			p.transform_point("satellite","degrees",lon,lat);
			break;
		case 7:
			p.transform_point("eckert iv","degrees",lon,lat);
			break;
		case 8:
			p.transform_point("equal earth","degrees",lon,lat);
			break;
		case 9:
			p.transform_point("epsg 3857","degrees",lon,lat);
			break;
	}
}

// ----------------------------------------------------------------------
void vtkOGSPointSource::SetDepth(double depth) {
	this->Center[2] = -depth*this->DepthScale;
	this->Modified();
}
void vtkOGSPointSource::GetDepth(double &depth) {
	depth = -this->Center[2] / this->DepthScale;
}