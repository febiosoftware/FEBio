// Interrupt.cpp: implementation of the Interrupt class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "Interrupt.h"
#include <signal.h>

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

bool Interruptable::m_bsig = false;

Interruptable::Interruptable()
{
	static bool binit = false;

	if (!binit) 
	{
		signal(SIGINT, Interruptable::handler);
		binit = true;
	}
}

Interruptable::~Interruptable()
{
	// TODO: restore original interruption handler
}

void Interruptable::handler(int sig)
{
	m_bsig = true;
	signal(SIGINT, Interruptable::handler);
}
