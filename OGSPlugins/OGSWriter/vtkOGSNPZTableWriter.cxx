/*=========================================================================

	Program:   OGSWriter
	Module:    vtkOGSNPZTableWriter.cxx

	Copyright (c) 2019 Arnau Miro, OGS
	All rights reserved.

		 This software is distributed WITHOUT ANY WARRANTY; without even
		 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
		 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkOGSNPZTableWriter.h"

#include "vtkCell.h"
#include "vtkPoints.h"
#include "vtkCallbackCommand.h"
#include "vtkCommand.h"
#include "vtkExecutive.h"
#include "vtkTable.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkFloatArray.h"
#include "vtkAbstractArray.h"
#include "vtkDoubleArray.h"
#include "vtkStringArray.h"

#include "vtkObjectFactory.h"

#include<vector>
#include<string>

vtkStandardNewMacro(vtkOGSNPZTableWriter);

//----------------------------------------------------------------------------

#include "macros.h"
#include "V3.h"
#include "field.h"
#include "vtkFields.h"
#include "vtkOperations.h"

#include "cnpy.hpp"

//----------------------------------------------------------------------------
vtkOGSNPZTableWriter::vtkOGSNPZTableWriter() : FileName(nullptr) {}

//----------------------------------------------------------------------------
vtkOGSNPZTableWriter::~vtkOGSNPZTableWriter() {
	this->SetFileName(nullptr);
}

//----------------------------------------------------------------------------
int vtkOGSNPZTableWriter::FillInputPortInformation( int vtkNotUsed(port), vtkInformation* info) {
	info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkTable");
	return 1;
}

//----------------------------------------------------------------------------
void vtkOGSNPZTableWriter::WriteData() {
	// Make sure we only export from one mesh
	if(this->GetNumberOfInputConnections(0) != 1)
		vtkErrorMacro("Exactly one input required.");

	// Recover the input
	vtkTable *input = vtkTable::SafeDownCast( this->GetExecutive()->GetInputData(0,0) );

	// Recover the number of rows and columns
	int nrows = input->GetNumberOfRows();
	int ncols = input->GetNumberOfColumns();

	// Loop the columns and save the data
	VTKARRAY *vtkcolumn; field::Field<FLDARRAY> column;
	for (int icol = 0; icol < ncols; ++icol) {
		vtkcolumn = VTKARRAY::SafeDownCast( input->GetColumn(icol) );
		column = VTK::createFieldfromVTK<VTKARRAY,FLDARRAY>(vtkcolumn);

		size_t Nn = (size_t)(column.get_n()), Nm = (size_t)(column.get_m());
		if (icol == 0) // Write
			cnpy::npz_save<FLDARRAY>(this->FileName,vtkcolumn->GetName(),column.data(),{Nn,Nm},"w");
		else // or append
			cnpy::npz_save<FLDARRAY>(this->FileName,vtkcolumn->GetName(),column.data(),{Nn,Nm},"a");
	}
}

//----------------------------------------------------------------------------
int vtkOGSNPZTableWriter::Write() {
	return Superclass::Write();
}