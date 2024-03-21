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



#include "stdafx.h"
#include "FEPlotData.h"

//-----------------------------------------------------------------------------
FEPlotData::FEPlotData(FEModel* fem) : FECoreBase(fem)
{
	m_ntype = PLT_FLOAT;
	m_sfmt = FMT_NODE;
	m_nregion = FE_REGION_NODE;

	m_arraySize = 0;
	m_szdom[0] = 0;
	m_szunit = nullptr;
}

//-----------------------------------------------------------------------------
FEPlotData::FEPlotData(FEModel* fem, Region_Type R, Var_Type t, Storage_Fmt s) : FECoreBase(fem)
{ 
	m_ntype = t; 
	m_sfmt = s; 
    m_nregion = R;

	m_arraySize = 0;
	m_szdom[0] = 0;

	m_szunit = nullptr;
}

//-----------------------------------------------------------------------------
void FEPlotData::SetArraySize(int n)
{
	m_arraySize = n;
}

//-----------------------------------------------------------------------------
int FEPlotData::GetArraysize() const
{
	return m_arraySize;
}

//-----------------------------------------------------------------------------
int FEPlotData::VarSize(Var_Type t)
{
	int ndata = 0;
	switch (DataType())
	{
	case PLT_FLOAT  : ndata =  1; break;
	case PLT_VEC3F  : ndata =  3; break;
	case PLT_MAT3FS : ndata =  6; break;
	case PLT_MAT3FD : ndata =  3; break;
    case PLT_TENS4FS: ndata = 21; break;
	case PLT_MAT3F  : ndata =  9; break;
	case PLT_ARRAY  : ndata = GetArraysize(); break;
	case PLT_ARRAY_VEC3F: ndata = GetArraysize()*3; break;
	}
	assert(ndata);
	return ndata;
}

//-----------------------------------------------------------------------------
void FEPlotData::SetDomainName(const char* szdom)
{
	strcpy(m_szdom, szdom); 
}

FEPlotFieldDescriptor::FEPlotFieldDescriptor(const std::string& typeString)
{
	m_valid = false;
	m_filterType = NO_FILTER;
	numFilter = -1;
	fieldName = typeString;

	// see if there is an alias defined
	size_t equalSign = fieldName.find('=');
	if (equalSign != string::npos)
	{
		alias = fieldName.substr(equalSign + 1, string::npos);
		fieldName.erase(equalSign, string::npos);
	}
	else alias = fieldName;

	// see if there is a filter
	size_t leftBracket = fieldName.find('[');
	if (leftBracket != string::npos)
	{
		// find the right bracket
		size_t rightBracket = fieldName.rfind(']');
		if (rightBracket == string::npos) return;

		string filter = fieldName.substr(leftBracket + 1, rightBracket - leftBracket - 1);
		fieldName.erase(leftBracket, string::npos);

		// see if the filter is a number or a string
		size_t leftQuote = filter.find('\'');
		if (leftQuote != string::npos)
		{
			size_t rightQuote = filter.rfind('\'');
			if ((rightQuote == string::npos) || (rightQuote == leftQuote)) return;

			m_filterType = STRING_FILTER;
			strFilter = filter.substr(leftQuote + 1, rightQuote - leftQuote - 1);
		}
		else
		{
			m_filterType = NUMBER_FILTER;
			numFilter = atoi(filter.c_str());
		}
	}

	m_valid = true;
}
