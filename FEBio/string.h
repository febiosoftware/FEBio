// string.h: interface for the string class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_STRING_H__D182FAC1_6408_469B_AF9F_C535E76CE45F__INCLUDED_)
#define AFX_STRING_H__D182FAC1_6408_469B_AF9F_C535E76CE45F__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <string.h>

class string  
{
public:
	string();
	string(const char* sz);
	string(const string& s);
	virtual ~string();

	string& operator = (const string& s);
	operator const char* () { return (const char*) m_sz; }

	string operator + (const string& s); // string concatentation
	string& operator +=(const string& s); // string concatentation

	int length() { return (m_sz?strlen(m_sz):0); }

	// comparison
	bool operator == (const string& s) { return (strcmp(m_sz, s.m_sz) == 0); }

protected:
	char*	m_sz;
};

#endif // !defined(AFX_STRING_H__D182FAC1_6408_469B_AF9F_C535E76CE45F__INCLUDED_)
