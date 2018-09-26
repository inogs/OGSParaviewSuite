/*=========================================================================

  Program:   ParaView
  Module:    vtkOGSVerticalProfilePlot.cxx

=========================================================================*/

#include "vtkOGSVerticalProfilePlot.h"

#include "vtkPythonRepresentation.h"

#include "vtkObjectFactory.h"

#include <sstream>
#include <string>

vtkStandardNewMacro(vtkOGSVerticalProfilePlot);

//----------------------------------------------------------------------------
vtkOGSVerticalProfilePlot::vtkOGSVerticalProfilePlot() {
	this->Script = NULL;
	strcpy(this->Params,"");
}

vtkOGSVerticalProfilePlot::~vtkOGSVerticalProfilePlot() {
	this->SetScript(NULL);
}

//----------------------------------------------------------------------------
void vtkOGSVerticalProfilePlot::Update() {
	// Enable all arrays for plotting
	this->Superclass::EnableAllAttributeArrays();

	// We append the python variables (in the form of name = value) to the script 
	std::ostringstream buf;
	buf << this->Params << "\n" << this->Script;
	// Update the Script property of the Superclass with the newly generated script
	this->Superclass::SetScript(buf.str().c_str());
	// Finally call the update method of the superclass
	this->Superclass::Update();
}

//----------------------------------------------------------------------------
void vtkOGSVerticalProfilePlot::SetParameterInternal(const char* raw_name, const char* raw_value) {
	const std::string name = raw_name ? raw_name : "";
	const std::string value = raw_value ? raw_value : "";

	if (name.empty()) {
		vtkErrorMacro(<< "cannot set parameter with empty name");
		return;
	}

	sprintf(this->Params,"%s%s = %s\n",this->Params,raw_name,raw_value);
	this->Modified();
}

void vtkOGSVerticalProfilePlot::SetParameter(const char* raw_name, const char* value) {
	std::ostringstream buf;

	buf << "r'" << value << "'";
	this->SetParameterInternal(raw_name, buf.str().c_str());
}

void vtkOGSVerticalProfilePlot::SetParameter(const char* raw_name, int value)
{
  std::ostringstream buf;
  buf << value;
  this->SetParameterInternal(raw_name, buf.str().c_str());
}

void vtkOGSVerticalProfilePlot::SetParameter(const char* raw_name, double value)
{
  std::ostringstream buf;
  buf << value;
  this->SetParameterInternal(raw_name, buf.str().c_str());
}

void vtkOGSVerticalProfilePlot::SetParameter(
  const char* raw_name, double value1, double value2, double value3)
{
  std::ostringstream buf;
  buf << "[" << value1 << ", " << value2 << ", " << value3 << "]";
  this->SetParameterInternal(raw_name, buf.str().c_str());
}

void vtkOGSVerticalProfilePlot::SetParameter(const char* raw_name, double value1, double value2)
{
  std::ostringstream buf;
  buf << "[" << value1 << ", " << value2 << "]";
  this->SetParameterInternal(raw_name, buf.str().c_str());
}