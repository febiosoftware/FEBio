// FileImport.cpp: implementation of the FEFileImport class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "FileImport.h"
#include <stdio.h>
#include <string.h>
#include <stdarg.h>

//-----------------------------------------------------------------------------
// FUNCTION : FEFileImport::FEFileImport
//
FEFileImport::FEFileImport()
{
	m_fp = 0;
	m_szerr[0] = 0;
}

//-----------------------------------------------------------------------------
// FUNCTION : FEFileImport::~FEFileImport
//
FEFileImport::~FEFileImport()
{
	// make sure to close the file
	Close();
}

//-----------------------------------------------------------------------------
// FUNCTION : FEFileImport::Open
//
bool FEFileImport::Open(const char* szfile, const char* szmode)
{
	m_fp = fopen(szfile, szmode);
	if (m_fp == 0) return false;

	strcpy(m_szfile, szfile);

	return true;
}

//-----------------------------------------------------------------------------
// FUNCTION : FEFileImport::Close
//
void FEFileImport::Close()
{
	if (m_fp)
	{
		fclose(m_fp);
		m_fp = 0;
	}
}

//-----------------------------------------------------------------------------
// FUNCTION : FEFileImport::GetErrorMessage
//
void FEFileImport::GetErrorMessage(char* szerr)
{
	strcpy(szerr, m_szerr);
}

//-----------------------------------------------------------------------------
// FUNCTION : FEFileImport::errf
//
bool FEFileImport::errf(const char* szerr, ...)
{
	// get a pointer to the argument list
	va_list	args;

	// copy to string
	va_start(args, szerr);
	vsprintf(m_szerr, szerr, args);
	va_end(args);

	// close the file
	Close();

	return false;
}
