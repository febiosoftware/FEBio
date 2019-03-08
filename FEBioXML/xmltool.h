#pragma once
#include "XMLReader.h"
#include <FECore/FEParameterList.h>
#include "febioxml_api.h"

//This namespace defines some helper functions that facilitate processing the FEBio xml formatted files.
namespace fexml
{
//---------------------------------------------------------------------------------------
// Reads the value of a parameter.
// if paramName is zero, the tag's name will be used as the parameter name.
bool FEBIOXML_API readParameter(XMLTag& tag, FEParameterList& paramList, const char* paramName = 0);

//---------------------------------------------------------------------------------------
// read a list of integers
void FEBIOXML_API readList(XMLTag& tag, vector<int>& l);
}
