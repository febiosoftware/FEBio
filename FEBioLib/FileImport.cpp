// FileImport.cpp: implementation of the FEFileImport class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "FileImport.h"
#include <stdio.h>
#include <string.h>
#include <stdarg.h>

//-----------------------------------------------------------------------------
//! class constructor
FEFileImport::FEFileImport()
{
	m_fp = 0;
	m_szerr[0] = 0;
}

//-----------------------------------------------------------------------------
// class destructor. Closes file on call.
FEFileImport::~FEFileImport()
{
	// make sure to close the file
	Close();
}

//-----------------------------------------------------------------------------
//! Open a file and store the file name and file pointer.
bool FEFileImport::Open(const char* szfile, const char* szmode)
{
	m_fp = fopen(szfile, szmode);
	if (m_fp == 0) return false;

	strcpy(m_szfile, szfile);

	return true;
}

//-----------------------------------------------------------------------------
//! Close the file
void FEFileImport::Close()
{
	if (m_fp)
	{
		fclose(m_fp);
		m_fp = 0;
	}
}

//-----------------------------------------------------------------------------
//! Get the error message. Errors message are stored when calling the errf function.
void FEFileImport::GetErrorMessage(char* szerr)
{
	strcpy(szerr, m_szerr);
}

//-----------------------------------------------------------------------------
//! Call this function to report an error message. The user can retrieve the 
//! error message with the GetErrorMessage member function.
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
