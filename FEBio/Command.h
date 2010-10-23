#pragma once

#include <string.h>

class Command
{
public:
	Command()
	{
		m_pszname = 0;
		m_pszdesc = 0;
	}
	virtual ~Command() {}

	void SetName(const char* szname)
	{
		int l = strlen(szname);
		assert(l);
		m_pszname = new char[l+1];
		strcpy(m_pszname, szname);
	}

	void SetDescription(const char* szdesc)
	{
		int l = strlen(szdesc);
		assert(l);
		m_pszdesc = new char[l+1];
		strcpy(m_pszdesc, szdesc);
	}

	const char* GetName() { return m_pszname; }
	const char* GetDescription() { return m_pszdesc; }

	virtual int run(int nargs, char** argv) = 0;

protected:
	char*	m_pszname;
	char*	m_pszdesc;
};
