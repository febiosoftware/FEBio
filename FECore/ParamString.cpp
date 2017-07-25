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
	string tmp;
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
bool ParamString::operator==(const string& s) const
{
	if (m_param.empty()) return false;
	return m_param[0]._name == s;
}

//-----------------------------------------------------------------------------
bool ParamString::operator!=(const string& s) const
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
