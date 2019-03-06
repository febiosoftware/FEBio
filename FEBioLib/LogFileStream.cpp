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
