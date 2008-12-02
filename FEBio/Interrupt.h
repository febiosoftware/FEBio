// Interrupt.h: interface for the Interrupt class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_INTERRUPT_H__AE370AC5_05F4_4290_B708_FFB4252F0AEB__INCLUDED_)
#define AFX_INTERRUPT_H__AE370AC5_05F4_4290_B708_FFB4252F0AEB__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

class Interruptable  
{
public:
	Interruptable();
	virtual ~Interruptable();

	static void handler(int sig);
	static bool	m_bsig;

	//! CTRL+C interruption handler
	void interrupt();
};

#endif // !defined(AFX_INTERRUPT_H__AE370AC5_05F4_4290_B708_FFB4252F0AEB__INCLUDED_)
