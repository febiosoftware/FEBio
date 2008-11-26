// console.h: interface for the Console class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_CONSOLE_H__CED42E3A_4BA6_44CE_8698_AB9C4328FB80__INCLUDED_)
#define AFX_CONSOLE_H__CED42E3A_4BA6_44CE_8698_AB9C4328FB80__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

class Console  
{
public:
	static void SetTitle(const char* sz, ...);

	static void Activate() { m_bActive = true; } 
	static void Deactivate() { m_bActive = false; }

protected:
	static	bool	m_bActive;
};

#endif // !defined(AFX_CONSOLE_H__CED42E3A_4BA6_44CE_8698_AB9C4328FB80__INCLUDED_)
