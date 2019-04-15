/*=========================================================================

  Program:   OGSPlotViews
  Module:    vtkOGSSpaghettiPlot.cxx

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkPythonRepresentation.h"

#include "vtkObjectFactory.h"

#include <sstream>
#include <string>
#include <iterator>

#include "vtkOGSSpaghettiPlot.h"

vtkStandardNewMacro(vtkOGSSpaghettiPlot);

//----------------------------------------------------------------------------
vtkOGSSpaghettiPlot::vtkOGSSpaghettiPlot() {
	this->Script    = NULL;
	this->Variables = NULL;
}

vtkOGSSpaghettiPlot::~vtkOGSSpaghettiPlot() {
	this->SetScript(NULL);
	this->SetVariables(NULL);
}

//----------------------------------------------------------------------------
void vtkOGSSpaghettiPlot::Update() {
	// Enable all arrays for plotting
	this->Superclass::EnableAllAttributeArrays();
	// First we iterate over the parameters map to create the variable = value pairs
	char aux1[2048]; strcpy(aux1,"");
	std::map<std::string,std::string>::iterator it;
	for (it = this->mapParam.begin(); it != mapParam.end(); it++)
		sprintf(aux1,"%s%s = %s\n",aux1,it->first.c_str(),it->second.c_str());
	// Second we create an auxiliary string where to store the plot variables
	char aux2[512];
	sprintf(aux2,"labels = '''%s'''",this->Variables);
	// We prepend both strings to the python script
	std::ostringstream buf;
	buf << aux1 << "\n" << aux2 << "\n" << this->Script;
	// Update the Script property of the Superclass with the newly generated script
	this->Superclass::SetScript(buf.str().c_str());
	// Finally call the update method of the superclass
	this->Superclass::Update();
}

//----------------------------------------------------------------------------
void vtkOGSSpaghettiPlot::SetParameterInternal(const char* raw_name, const char* raw_value) {
	const std::string name = raw_name ? raw_name : "";
	const std::string value = raw_value ? raw_value : "";

	if (name.empty()) {
		vtkErrorMacro(<< "cannot set parameter with empty name");
		return;
	}

	// Check if we can insert the parameter in the map
	if (!this->mapParam.insert(std::make_pair(raw_name,raw_value)).second) {
		// We couldn't insert the key, which means it already exists
		// we shall update its value
		mapParam[raw_name] = raw_value;
	}

	this->Modified();
}

void vtkOGSSpaghettiPlot::SetParameter(const char* raw_name, const char* value) {
	std::ostringstream buf;

	buf << "r'" << value << "'";
	this->SetParameterInternal(raw_name, buf.str().c_str());
}

void vtkOGSSpaghettiPlot::SetParameter(const char* raw_name, int value)
{
  std::ostringstream buf;
  buf << value;
  this->SetParameterInternal(raw_name, buf.str().c_str());
}

void vtkOGSSpaghettiPlot::SetParameter(const char* raw_name, double value)
{
  std::ostringstream buf;
  buf << value;
  this->SetParameterInternal(raw_name, buf.str().c_str());
}

void vtkOGSSpaghettiPlot::SetParameter(
  const char* raw_name, double value1, double value2, double value3)
{
  std::ostringstream buf;
  buf << "[" << value1 << ", " << value2 << ", " << value3 << "]";
  this->SetParameterInternal(raw_name, buf.str().c_str());
}

void vtkOGSSpaghettiPlot::SetParameter(const char* raw_name, double value1, double value2)
{
  std::ostringstream buf;
  buf << "[" << value1 << ", " << value2 << "]";
  this->SetParameterInternal(raw_name, buf.str().c_str());
}