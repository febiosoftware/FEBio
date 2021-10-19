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
#include "ParamString.h"

//=============================================================================
void cstrncpy(char* pd, char* ps, int l)
{
	for (int i = 0; i <= l; ++i) pd[i] = ps[i];
}

ParamString::ParamString(const char* szparam)
{
	ParamRef p;
	std::string tmp;
	const char* sz = szparam;
	int in = 0;
	while (*sz)
	{
		if (*sz == '.')
		{
			// add parameter
			m_param.push_back(p);

			// reset
			p.clear();
			tmp.clear();
			in = 0;
		}
		else if (*sz == '[')
		{
			in = 1;
		}
		else if (*sz == ']')
		{
			if (in == 1)
			{
				const char* s = tmp.c_str();
				p._index = atoi(s);
			}
		}
		else if (*sz == '(')
		{
			in = 2;
		}
		else if (*sz == ')')
		{
			if (in == 2)
			{
				const char* s = tmp.c_str();
				p._id = atoi(s);
			}
		}
		else if (*sz == '\'')
		{
			if (in == 2)
			{
				in = 3;
			}
		}
		else if (*sz == '"')
		{
			if (in == 2)
			{
				in = 3;
			}
		}
		else 
		{
			if      (in == 0) p._name.push_back(*sz);
			else if (in == 1) tmp.push_back(*sz);
			else if (in == 2) tmp.push_back(*sz);
			else if (in == 3) p._idName.push_back(*sz);
		}

		// next char
		sz++;
	}

	// don't forget to add the last ref
	m_param.push_back(p);
}

//-----------------------------------------------------------------------------
ParamString::ParamString(const ParamString& p)
{
	m_param = p.m_param;
}

//-----------------------------------------------------------------------------
void ParamString::operator=(const ParamString& p)
{
	m_param = p.m_param;
}

//-----------------------------------------------------------------------------
ParamString::~ParamString()
{
}

//-----------------------------------------------------------------------------
//! number of refs in string
int ParamString::count() const
{
	return (int) m_param.size();
}

//-----------------------------------------------------------------------------
ParamString ParamString::next() const
{
	ParamString p;
	if (m_param.size() > 1)
	{
		p.m_param.insert(p.m_param.end(), m_param.begin() + 1, m_param.end());
	}
	return p;
}

//-----------------------------------------------------------------------------
bool ParamString::isValid() const { return (m_param.empty() == false); }

//-----------------------------------------------------------------------------
void ParamString::clear() { m_param.clear(); }

//-----------------------------------------------------------------------------
ParamString ParamString::last() const
{
	ParamString p;
	if (m_param.size() > 1)
	{
		p.m_param.push_back(m_param[m_param.size()-1]);
	}
	return p;
}

//-----------------------------------------------------------------------------
bool ParamString::operator==(const std::string& s) const
{
	if (m_param.empty()) return false;
	return m_param[0]._name == s;
}

//-----------------------------------------------------------------------------
bool ParamString::operator!=(const std::string& s) const
{
	if (m_param.empty()) return false;
	return m_param[0]._name != s;
}

//-----------------------------------------------------------------------------
const char* ParamString::c_str() const
{
	if (m_param.empty()) return 0;
	return m_param[0]._name.c_str();
}

//-----------------------------------------------------------------------------
std::string ParamString::string() const
{
	if (m_param.empty()) return 0;
	return m_param[0]._name;
}

//-----------------------------------------------------------------------------
//! get the index (-1 if index not a number)
int ParamString::Index() const
{
	if (m_param.empty()) return -1;
	return m_param[0]._index;
}

//-----------------------------------------------------------------------------
//! get the ID (-1 if index not a number)
int ParamString::ID() const
{
	if (m_param.empty()) return -1;
	return m_param[0]._id;
}

//-----------------------------------------------------------------------------
//! get the index name (null if not defined)
const char* ParamString::IDString() const
{
	if (m_param.empty()) return 0;
	if (m_param[0]._idName.empty()) return 0;
	return m_param[0]._idName.c_str();
}
