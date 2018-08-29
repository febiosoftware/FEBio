#pragma once
#include "XMLReader.h"
#include <FECore/FEParameterList.h>

//This namespace defines some helper functions that facilitate processing the FEBio xml formatted files.
namespace fexml
{
//---------------------------------------------------------------------------------------
// Reads the value of a parameter.
// if paramName is zero, the tag's name will be used as the parameter name.
bool readParameter(XMLTag& tag, FEParameterList& paramList, const char* paramName = 0);

//---------------------------------------------------------------------------------------
// read a list of integers
void readList(XMLTag& tag, vector<int>& l);
}
