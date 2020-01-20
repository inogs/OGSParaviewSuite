/*=========================================================================

	Program:   OGSUtils
	Module:    vtkOGSDepthLineSource.cxx

	Copyright (c) 2018 Arnau Miro, OGS
	All rights reserved.

		 This software is distributed WITHOUT ANY WARRANTY; without even
		 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
		 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkOGSDepthLineSource.h"

#include "vtkCellArray.h"
#include "vtkFloatArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"

#include "projection.h"

#include <cmath>

vtkStandardNewMacro(vtkOGSDepthLineSource);
vtkCxxSetObjectMacro(vtkOGSDepthLineSource,Points,vtkPoints);

// ----------------------------------------------------------------------
vtkOGSDepthLineSource::vtkOGSDepthLineSource(int res) {

	this->Point1[0] = -.5;
	this->Point1[1] =  .0;
	this->Point1[2] =  .0;

	this->Point2[0] =  .5;
	this->Point2[1] =  .0;
	this->Point2[2] =  .0;

	this->DepthScale = 1000.;
	this->Projection = 0;

	this->Points = nullptr;

	this->Resolution = ( res < 1 ? 1 : res );
	this->OutputPointsPrecision = SINGLE_PRECISION;

	this->SetNumberOfInputPorts( 0 );
}

// ----------------------------------------------------------------------
vtkOGSDepthLineSource::~vtkOGSDepthLineSource() {

	this->SetPoints( nullptr );
}

// ----------------------------------------------------------------------
int vtkOGSDepthLineSource::RequestInformation(vtkInformation *vtkNotUsed(request),
	vtkInformationVector **vtkNotUsed(inputVector), vtkInformationVector *outputVector) {

	// get the info object
	vtkInformation *outInfo = outputVector->GetInformationObject(0);
	outInfo->Set(CAN_HANDLE_PIECE_REQUEST(),1);
	return 1;
}

// ----------------------------------------------------------------------
int vtkOGSDepthLineSource::RequestData(vtkInformation *vtkNotUsed(request),
	vtkInformationVector **vtkNotUsed(inputVector), vtkInformationVector *outputVector) {

	// Reject meaningless parameterizations
	vtkIdType nSegments = this->Points ? this->Points->GetNumberOfPoints() - 1 : 1;
	if ( nSegments < 1 ) {
		vtkWarningMacro( <<"Cannot define a broken line with given input.");
		return 0;
	}

	// get the info object
	vtkInformation *outInfo = outputVector->GetInformationObject(0);

	// get the output
	vtkPolyData *output = vtkPolyData::SafeDownCast(
		outInfo->Get(vtkDataObject::DATA_OBJECT()));

	if (outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER()) > 0)
		return 1;

	// Create and allocate lines
	vtkIdType numLines = nSegments * this->Resolution;
	vtkCellArray *newLines = vtkCellArray::New();
	newLines->Allocate( newLines->EstimateSize( numLines, 2 ) );

	// Create and allocate points
	vtkIdType numPts = numLines + 1;
	vtkPoints *newPoints = vtkPoints::New();

	// Set the desired precision for the points in the output.
	if(this->OutputPointsPrecision == vtkAlgorithm::DOUBLE_PRECISION)
		newPoints->SetDataType(VTK_DOUBLE);
	else
		newPoints->SetDataType(VTK_FLOAT);

	newPoints->Allocate( numPts );

	// Create and allocate texture coordinates
	vtkFloatArray *newTCoords = vtkFloatArray::New();
	newTCoords->SetNumberOfComponents( 2 );
	newTCoords->Allocate( 2 * numPts );
	newTCoords->SetName( "Texture Coordinates" );

	// Allocate convenience storage
	double x[3], tc[3], v[3];

	// Generate points and texture coordinates
	if ( this->Points ) {

		// Create storage for segment endpoints
		double point1[3];
		double point2[3];

		// Point index offset for fast insertion
		vtkIdType offset = 0;

		// Iterate over segments
		for ( vtkIdType s = 0; s < nSegments; ++ s ) {

			// Get coordinates of endpoints
			this->Points->GetPoint( s, point1 );
			this->Points->GetPoint( s + 1, point2 );

			// Calculate segment vector
			for ( int i = 0; i < 3; ++ i )
			{
				v[i] = point2[i] - point1[i];
			}

			// Generate points along segment
			tc[1] = 0.;
			tc[2] = 0.;
			for ( vtkIdType i = 0; i < this->Resolution; ++ i, ++ offset ) {
				tc[0] = static_cast<double>( i ) / this->Resolution;
				for ( int j = 0; j < 3; ++ j )
					x[j] = point1[j] + tc[0] * v[j];
				newPoints->InsertPoint( offset, x );
				newTCoords->InsertTuple( offset, tc );
			}
		} // s

		// Generate last endpoint
		newPoints->InsertPoint( numLines, point2 );
		tc[0] = 1.;
		newTCoords->InsertTuple( numLines, tc );

	} else {
		// Calculate segment vector
		for ( int i = 0; i < 3; ++ i )
			v[i] = this->Point2[i] - this->Point1[i];

		// Generate points along segment
		tc[1] = 0.;
		tc[2] = 0.;
		for ( vtkIdType i = 0; i < numPts; ++ i ) {
			tc[0] = static_cast<double>( i ) / this->Resolution;
			for ( int j = 0; j < 3; ++ j )
				x[j] = this->Point1[j] + tc[0] * v[j];
			newPoints->InsertPoint( i, x );
			newTCoords->InsertTuple( i, tc );
		}
	} // else

	//  Generate lines
	newLines->InsertNextCell( numPts );
	for ( vtkIdType i = 0; i < numPts; ++ i )
		newLines->InsertCellPoint( i );

	// Update ourselves and release memory
	output->SetPoints( newPoints );
	newPoints->Delete();

	output->GetPointData()->SetTCoords( newTCoords );
	newTCoords->Delete();

	output->SetLines( newLines );
	newLines->Delete();

	return 1;
}

// ----------------------------------------------------------------------
void vtkOGSDepthLineSource::SetLonLat(double lon, double lat) {
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

	// Set points
	this->Point1[0] = this->Point2[0] = lon;
	this->Point1[1] = this->Point2[1] = lat;
	this->Modified();
}
void vtkOGSDepthLineSource::GetLonLat(double &lon, double &lat) {
	// Conversion to lon, lat
	PROJ::Projection p;
	lon = this->Point1[0]; lat = this->Point1[2];
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
void vtkOGSDepthLineSource::SetDepthScale(double d) {
	// Recover depth range
	double d1,d2;
	this->GetDepthRange(d1,d2);
	// Set new depth scale
	this->DepthScale = d;
	// Set new depth range
	this->SetDepthRange(d1,d2);
	this->Modified();
}

// ----------------------------------------------------------------------
void vtkOGSDepthLineSource::SetDepthRange(double d1, double d2) {
	this->Point1[2] = d1*this->DepthScale;
	this->Point2[2] = -d2*this->DepthScale;
	this->Modified();
}
void vtkOGSDepthLineSource::GetDepthRange(double &d1, double &d2) {
	d1 = this->Point1[2]/this->DepthScale;
	d2 = -this->Point2[2]/this->DepthScale;
}

// ----------------------------------------------------------------------
void vtkOGSDepthLineSource::SetPoint1(double p0, double p1, double p2) {
	this->Point1[0] = this->Point2[0] = p0;
	this->Point1[1] = this->Point2[1] = p1;
	this->Point1[2] = p2;
	this->Modified();
}
void vtkOGSDepthLineSource::SetPoint1(double p[3]) {
	SetPoint1(p[0],p[1],p[2]);
}
void vtkOGSDepthLineSource::SetPoint1(float point1f[3]) {
	double point1d[3];
	point1d[0] = point1f[0];
	point1d[1] = point1f[1];
	point1d[2] = point1f[2];
	SetPoint1(point1d);
}
void vtkOGSDepthLineSource::GetPoint1(double &p0, double &p1, double &p2) {
	p0 = this->Point1[0];
	p1 = this->Point1[1];
	p2 = this->Point1[2];
}
void vtkOGSDepthLineSource::GetPoint1(double p[3]) {
	GetPoint1(p[0],p[1],p[2]);
}

// ----------------------------------------------------------------------
void vtkOGSDepthLineSource::SetPoint2(double p0, double p1, double p2) {
	this->Point2[0] = this->Point1[0];
	this->Point2[1] = this->Point1[1];
	this->Point2[2] = p2;
	this->Modified();
}
void vtkOGSDepthLineSource::SetPoint2(double point[3]) {
	SetPoint1(point[0],point[1],point[2]);
}
void vtkOGSDepthLineSource::SetPoint2(float point2f[3]) {
	double point2d[3];
	point2d[0] = point2f[0];
	point2d[1] = point2f[1];
	point2d[2] = point2f[2];
	SetPoint2(point2d);
}
void vtkOGSDepthLineSource::GetPoint2(double &p0, double &p1, double &p2) {
	p0 = this->Point2[0];
	p1 = this->Point2[1];
	p2 = this->Point2[2];
}
void vtkOGSDepthLineSource::GetPoint2(double p[3]) {
	GetPoint2(p[0],p[1],p[2]);
}
