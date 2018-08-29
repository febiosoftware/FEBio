// Interrupt.cpp: implementation of the Interrupt class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "Interrupt.h"
#include "CommandManager.h"
#include "console.h"
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
//! \todo Restore original intteruption handler
Interruption::~Interruption()
{
	
}

void Interruption::handler(int sig)
{
	m_bsig = true;
	signal(SIGINT, Interruption::handler);
}

void Interruption::interrupt()
{
	// get a pointer to the console window
	Console* pShell = Console::GetHandle();

	// get a pointer to the command manager
	CommandManager* pCM = CommandManager::GetInstance();

	int nargs;
	char* argv[32];

	// enter command loop
	while (1)
	{
		// get a command from the shell
		pShell->GetCommand(nargs, argv);
		if (nargs > 0)
		{
			// find the command that has this name
			Command* pcmd = pCM->Find(argv[0]);
			if (pcmd)
			{
				int nret = pcmd->run(nargs, argv);
				if (nret == 1) break;
			}
			else
			{
				printf("Unknown command: %s\n", argv[0]);
			}
		}
	}
}
