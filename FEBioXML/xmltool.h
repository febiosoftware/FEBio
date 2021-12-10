/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
the City of New York, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/



#pragma once
#include <XML/XMLReader.h>
#include <FECore/FECoreBase.h>
#include <FECore/ClassDescriptor.h>
#include "febioxml_api.h"

//This namespace defines some helper functions that facilitate processing the FEBio xml formatted files.
namespace fexml
{
//---------------------------------------------------------------------------------------
// Reads the value of a parameter.
// if paramName is zero, the tag's name will be used as the parameter name.
bool FEBIOXML_API readParameter(XMLTag& tag, FEParameterList& paramList, const char* paramName = 0);

//---------------------------------------------------------------------------------------
bool FEBIOXML_API readParameter(XMLTag& tag, FECoreBase* pc);

//---------------------------------------------------------------------------------------
// reads the parameters and properties of a FECore class
bool FEBIOXML_API readParameterList(XMLTag& tag, FECoreBase* pc);

//---------------------------------------------------------------------------------------
// read a list of integers
void FEBIOXML_API readList(XMLTag& tag, vector<int>& l);

//---------------------------------------------------------------------------------------
// create a class descriptor from the current tag
FEBIOXML_API FEClassDescriptor* readParameterList(XMLTag& tag);

}
