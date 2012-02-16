// Document.cpp: implementation of the CDocument class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "Document.h"

//-----------------------------------------------------------------------------
CDocument::CDocument()
{
}

//-----------------------------------------------------------------------------
CDocument::~CDocument()
{
}

//-----------------------------------------------------------------------------
bool CDocument::OpenFile(const char* szfile)
{
	return m_fem.Input(szfile);
}
