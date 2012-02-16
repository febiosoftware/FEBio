// Document.h: interface for the CDocument class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_DOCUMENT_H__E699CE33_76BC_46FB_8CFC_4FA83D106B4C__INCLUDED_)
#define AFX_DOCUMENT_H__E699CE33_76BC_46FB_8CFC_4FA83D106B4C__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "FEM.h"

//-----------------------------------------------------------------------------
// Document class - stores all data
//
class CDocument  
{
public:
	// constructor/destructor
	CDocument();
	virtual ~CDocument();

protected:
	FEM		m_fem;
};

#endif // !defined(AFX_DOCUMENT_H__E699CE33_76BC_46FB_8CFC_4FA83D106B4C__INCLUDED_)
