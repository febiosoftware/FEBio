// string.cpp: implementation of the string class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "string.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

string::string()
{
	m_sz = 0;
}

string::string(const char* sz)
{
	if (sz)
	{
		int l = strlen(sz);
		m_sz = new char[l+1];
		strcpy(m_sz, sz);
	}
	else m_sz = 0;
}

string::string(const string& s)
{
	if (s.m_sz)
	{
		int l = strlen(s.m_sz);
		m_sz = new char[l+1];
		strcpy(m_sz, s.m_sz);
	}
	else m_sz = 0;
}

string& string::operator =(const string& s)
{
	if (m_sz) delete [] m_sz;
	m_sz = 0;

	if (s.m_sz)
	{
		int l = strlen(s.m_sz);
		m_sz = new char[l+1];
		strcpy(m_sz, s.m_sz);
	}
	
	return (*this);
}

string::~string()
{
	if (m_sz) delete [] m_sz;
}


string string::operator +(const string& s)
{
	int l0 = (m_sz?strlen(m_sz):0);
	int l1 = (s.m_sz?strlen(s.m_sz):0);

	int l = l0+l1;

	string sr;
	if (l > 0)
	{
		sr.m_sz = new char[l+1];
		if (m_sz) strcpy(sr.m_sz, m_sz);
		if (s.m_sz) strcpy(sr.m_sz+l0, s.m_sz);
	}

	return sr;
}

string& string::operator +=(const string& s)
{
	int l0 = (m_sz?strlen(m_sz):0);
	int l1 = (s.m_sz?strlen(s.m_sz):0);

	int l = l0+l1;

	if (l > 0)
	{
		char* sz = new char[l+1];
		if (m_sz) strcpy(sz, m_sz);
		if (s.m_sz) strcpy(sz+l0, s.m_sz);
		if (m_sz) delete [] m_sz;
		m_sz = sz;
	}

	return (*this);
}
