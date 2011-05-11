// FileImport.h: interface for the FileImport class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FILEIMPORT_H__AC15F5F8_E069_4640_B3FD_077984EEA78F__INCLUDED_)
#define AFX_FILEIMPORT_H__AC15F5F8_E069_4640_B3FD_077984EEA78F__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <stdio.h>
#include "fem.h"

//-----------------------------------------------------------------------------
//! Base class for file import classes

class FEFileImport  
{
public:
	FEFileImport();
	virtual ~FEFileImport();

	// This function will be overloaded in the derived classes
	virtual bool Load(FEM& fem, const char* szfile) = 0;

	void GetErrorMessage(char* szerr);

protected:
	bool Open(const char* szfile, const char* szmode);
	void Close();

	bool errf(const char* szerr, ...);

protected:
	FILE*	m_fp;
	char	m_szfile[256];	//!< file name

	char	m_szerr[256];	//!< error message
};


#endif // !defined(AFX_FILEIMPORT_H__AC15F5F8_E069_4640_B3FD_077984EEA78F__INCLUDED_)
