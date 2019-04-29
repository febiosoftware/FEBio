/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, The Trustees of Columbia University in 
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



#include "stdafx.h"
#include "xmltool.h"


//-----------------------------------------------------------------------------
//! This function parses a parameter list
bool fexml::readParameter(XMLTag& tag, FEParameterList& paramList, const char* paramName)
{
	// see if we can find this parameter
	FEParam* pp = paramList.FindFromName((paramName == 0 ? tag.Name() : paramName));
	if (pp == 0) return false;
	
	if (pp->dim() == 1)
	{
		switch (pp->type())
		{
		case FE_PARAM_DOUBLE  : { double d; tag.value(d); pp->value<double  >() = d; } break;
		case FE_PARAM_INT     : { int n; tag.value(n); pp->value<int>() = n; } break;
		case FE_PARAM_BOOL    : { bool b; tag.value(b); pp->value<bool>() = b; } break;
		default:
			assert(false);
			return false;
		}
	}
	else
	{
		switch (pp->type())
		{
		case FE_PARAM_INT   : { vector<int> d(pp->dim()); tag.value(&d[0], pp->dim()); } break;
		case FE_PARAM_DOUBLE: { vector<double> d(pp->dim()); tag.value(&d[0], pp->dim()); } break;
        default: break;
		}
	}

	return true;
}

//-----------------------------------------------------------------------------
void fexml::readList(XMLTag& tag, vector<int>& l)
{
	// make sure the list is empty
	l.clear();

	// get a pointer to the value
	const char* sz = tag.szvalue();

	// parse the string
	const char* ch;
	do
	{
		int n0, n1, nn;
		int nread = sscanf(sz, "%d:%d:%d", &n0, &n1, &nn);
		switch (nread)
		{
		case 1:
			n1 = n0;
			nn = 1;
			break;
		case 2:
			nn = 1;
			break;
		}

		for (int i = n0; i <= n1; i += nn) l.push_back(i);

		ch = strchr(sz, ',');
		if (ch) sz = ch + 1;
	}
	while (ch != 0);
}
