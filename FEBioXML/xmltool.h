#pragma once
#include "XMLReader.h"
#include <FECore/FEParameterList.h>

//---------------------------------------------------------------------------------------
// Reads the value of a parameter.
// if paramName is zero, the tag's name will be used as the parameter name.
bool xmlReadParameter(XMLTag& tag, FEParameterList& paramList, const char* paramName = 0);
