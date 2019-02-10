// Interrupt.cpp: implementation of the Interrupt class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "Interrupt.h"
#include "FEBioApp.h"
#include <signal.h>

bool Interruption::m_bsig = false;

Interruption::Interruption()
{
	static bool binit = false;

	if (!binit) 
	{
		signal(SIGINT, Interruption::handler);
		binit = true;
	}
}

//-----------------------------------------------------------------------------
//! Destructor
// TODO: Restore original intteruption handler
Interruption::~Interruption()
{
	
}

void Interruption::handler(int sig)
{
	m_bsig = true;
	signal(SIGINT, Interruption::handler);
}
