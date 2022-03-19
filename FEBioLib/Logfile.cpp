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
#include "Logfile.h"
#include <stdarg.h>
#include <cstring>

//-----------------------------------------------------------------------------
// constructor for the Logfile class
Logfile::Logfile()
{
	m_fp = 0;
	m_ps = 0;

	m_mode = LOG_FILE_AND_SCREEN;
}

//-----------------------------------------------------------------------------
// destructor for the Logfile class
//
Logfile::~Logfile()
{
	close();
}

//-----------------------------------------------------------------------------
// open a file
//
bool Logfile::open(const char* szfile)
{
	if (m_fp == 0) m_fp = new LogFileStream;
	return m_fp->open(szfile);
}

//-----------------------------------------------------------------------------
//  opens a file and prepares for appending
//
bool Logfile::append(const char* szfile)
{
	// store a copy of the filename
	if (m_fp == 0) m_fp = new LogFileStream;
	return m_fp->append(szfile);
}

//-----------------------------------------------------------------------------
//! flush the logfile
void Logfile::flush()
{
	if (m_fp) m_fp->flush(); 
	if (m_ps) m_ps->flush();
}

//-----------------------------------------------------------------------------
//! close the logfile
void Logfile::close()
{ 
	if (m_fp)
	{
		m_fp->close(); 
		delete m_fp;
		m_fp = 0;
	}
}

//-----------------------------------------------------------------------------
// This function works like all other printf functions
// with the exception that everything that is output to the file
// is (optionally) also output to the screen.
//
void Logfile::printf(const char* sz, ...)
{
	// get a pointer to the argument list
	va_list	args;

	// make the message
	char sztxt[2048] = {0};
	va_start(args, sz);
	vsprintf(sztxt, sz, args);
	va_end(args);
	
	// print to file
	if (m_fp && (m_mode & LOG_FILE)) m_fp->print(sztxt);

	// print to screen
	if (m_ps && (m_mode & LOG_SCREEN)) m_ps->print(sztxt);
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
	char sztxt[1024] = {0};
	va_start(args, sz);
	vsprintf(sztxt, sz, args);
	va_end(args);

	// print the box
	char szmsg[1024] = {0};
	char* ch = szmsg;
	sprintf(szmsg,"\n *************************************************************************\n"); ch += strlen(ch);
	// print the title
	if (sztitle)
	{
		int l = (int)strlen(sztitle);
		char left[60] = {0};
		char right[60] = {0};
		strncpy(left, sztitle, l/2);
		strncpy(right, sztitle+l/2, l - l/2);
		sprintf(ch," * %33s", left); ch += strlen(ch);
		sprintf(ch,"%-36s *\n", right); ch += strlen(ch);
//		sprintf(ch," *%71s*\n", ""); ch += strlen(ch);
	}

	// print the message
	char* ct = sztxt, *cn;
	char tmp;
	do
	{
		cn = strchr(ct,'\n');
		if (cn) *cn = 0;
		int l = (int)strlen(ct);
		bool wrap = false;
		int n = 69;
		if (l > 69) {
			while ((n > 0) && (ct[n] != ' ')) n--;
			if (n == 0) n = 69;
			tmp = ct[n];
			ct[n] = 0;
			wrap = true;
			if (cn) *cn = '\n';
			cn = ct + n; 
		}
		sprintf(ch," * %-69s *\n", ct); ch += strlen(ch);
		if (wrap) { ct[n] = tmp; }
		if (cn) ct = cn+1;
	}
	while (cn);
//	sprintf(ch," *                                                                       *\n"); ch += strlen(ch);
	sprintf(ch," *************************************************************************\n");

	// print the message
	printf(szmsg);
}


//-----------------------------------------------------------------------------
// Sets the Logfile mode. 
//
Logfile::MODE Logfile::SetMode(Logfile::MODE mode)
{
	MODE old = m_mode;
	m_mode = mode;
	return old;
}

//! get the loggin mode
Logfile::MODE Logfile::GetMode() 
{ 
	return m_mode; 
}
