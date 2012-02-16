// Wnd.h: interface for the CWnd class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_WND_H__793D79A3_EBE4_4660_8EE3_0016B7467520__INCLUDED_)
#define AFX_WND_H__793D79A3_EBE4_4660_8EE3_0016B7467520__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <Flx_Wnd.h>
#include "Document.h"

class CWnd : public Flx_Wnd
{
public:
	CWnd(int w, int h, const char* sztitle, CDocument* pdoc);
	virtual ~CWnd();

	CDocument* GetDocument() { return m_pDoc; }

	bool OpenFile(const char* szfile);

protected:
	CDocument*	m_pDoc;
};

#endif // !defined(AFX_WND_H__793D79A3_EBE4_4660_8EE3_0016B7467520__INCLUDED_)
