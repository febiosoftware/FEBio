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
#include "LogFileStream.h"

//-----------------------------------------------------------------------------
LogFileStream::LogFileStream()
{
	m_fp = 0;
}

//-----------------------------------------------------------------------------
LogFileStream::~LogFileStream()
{
	close();
}

//-----------------------------------------------------------------------------
void LogFileStream::close()
{
	if (m_fp) fclose(m_fp);
	m_fp = 0;
}

//-----------------------------------------------------------------------------
void LogFileStream::flush()
{
	if (m_fp) fflush(m_fp);
}

//-----------------------------------------------------------------------------
bool LogFileStream::open(const char* szfile)
{
	m_fileName = szfile;
	if (m_fp) close();
	m_fp = fopen(szfile, "wt");
	return (m_fp != NULL);
}

//-----------------------------------------------------------------------------
bool LogFileStream::append(const char* szfile)
{
	// make sure we don't have a log file already open
	if (m_fp)
	{
		fseek(m_fp, SEEK_END, 0);
		return true;
	}

	// create the log file
	m_fileName = szfile;
	m_fp = fopen(szfile, "a+t");

	return (m_fp != NULL);
}

//-----------------------------------------------------------------------------
void LogFileStream::print(const char* sztxt)
{
	if (m_fp) fprintf(m_fp, "%s", sztxt);
}
