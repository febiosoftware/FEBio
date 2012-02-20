// Logfile.cpp: implementation of the Logfile class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "Logfile.h"
#include <stdarg.h>
#include <string.h>

//-----------------------------------------------------------------------------
// The one-and-only logfile
Logfile& clog = *Logfile::GetInstance();

//-----------------------------------------------------------------------------
Logfile* Logfile::m_plog = 0;

//-----------------------------------------------------------------------------

Logfile* Logfile::GetInstance()
{
	if (m_plog == 0) m_plog = new Logfile();
	return m_plog;
}

//-----------------------------------------------------------------------------
// FUNCTION: Logfile::Logfile
// constructor for the Logfile class
//
Logfile::Logfile()
{
	m_fp = 0;
	m_szfile[0] = 0;

	m_mode = FILE_AND_SCREEN;
}

//-----------------------------------------------------------------------------
// FUNCTION: Logfile::~Logfile
// destructor for the Logfile class
//
Logfile::~Logfile()
{
	if (m_fp) fclose(m_fp);
	m_plog = 0;
}

//-----------------------------------------------------------------------------
// FUNCTION: Logfile::open
// open a file
//
bool Logfile::open(const char* szfile)
{
	// create the log file
	m_fp = fopen(szfile, "wt");

	strcpy(m_szfile, szfile);

	return (m_fp != 0);
}

//-----------------------------------------------------------------------------
// FUNCTION: Logfile::append
//  opens a file and prepares for appending
//

bool Logfile::append(const char* szfile)
{
	// make sure we don't have a log file already open
	if (m_fp) return false;

	// create the log file
	m_fp = fopen(szfile, "a+t");

	// store a copy of the filename
	strcpy(m_szfile, szfile);

	return (m_fp != 0);
}


//-----------------------------------------------------------------------------
// FUNCTION: Logfile::print
// This function works like all other printf functions
// with the exception that everything that is output to the file
// is (optionally) also output to the screen.
//
void Logfile::printf(const char* sz, ...)
{
	static char szmsg[1024] = {0};

	// get a pointer to the argument list
	va_list	args;

	// make the message
	static char sztxt[1024] = {0};
	va_start(args, sz);
	vsprintf(sztxt, sz, args);
	va_end(args);
	
	// print to file
	if (m_fp && (m_mode & FILE_ONLY)) fprintf(m_fp, sztxt);

	// print to screen
	if (m_ps && (m_mode & SCREEN_ONLY)) m_ps->print(sztxt);
}

//-----------------------------------------------------------------------------
// FUNCTION: Logfile::printbox
// This function prints a message insided a box
//
void Logfile::printbox(const char* sztitle, const char* sz, ...)
{
	// get a pointer to the argument list
	va_list	args;

	// make the message
	static char sztxt[1024] = {0};
	va_start(args, sz);
	vsprintf(sztxt, sz, args);
	va_end(args);

	// print the box
	char szmsg[1024] = {0};
	char* ch = szmsg;
	sprintf(szmsg," *************************************************************************\n"); ch += strlen(ch);
	// print the title
	if (sztitle)
	{
		int l = strlen(sztitle);
		char left[60] = {0};
		char right[60] = {0};
		strncpy(left, sztitle, l/2);
		strncpy(right, sztitle+l/2, l - l/2);
		sprintf(ch," * %33s", left); ch += strlen(ch);
		sprintf(ch,"%-36s *\n", right); ch += strlen(ch);
		sprintf(ch," *%71s*\n", ""); ch += strlen(ch);
	}

	// print the message
	char* ct = sztxt, *cn;
	do
	{
		cn = strchr(ct,'\n');
		if (cn) *cn = 0;
		sprintf(ch," * %-69s *\n", ct); ch += strlen(ch);
		if (cn) ct = cn+1;
	}
	while (cn);
	sprintf(ch," *************************************************************************\n");

	// print the message
	printf(szmsg);
}


//-----------------------------------------------------------------------------
// FUNCTION: Logfile::SetMode
// Sets the Logfile mode. 
//
Logfile::MODE Logfile::SetMode(Logfile::MODE mode)
{
	MODE old = m_mode;
	m_mode = mode;
	return old;
}
