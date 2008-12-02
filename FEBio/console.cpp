// console.cpp: implementation of the Console class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "console.h"
#include <stdarg.h>

#ifdef WIN32
	#include <windows.h>
#endif

//--------------------------------------------------------------------
// pointer to the one and only console object. This pointer
// will be initialized during the first call to GetHandle.
Console* Console::m_pShell = 0;

//--------------------------------------------------------------------
//! This class returns the pointer to the console object. On the first
//! call, the pointer is allocated.
Console* Console::GetHandle()
{
	if (m_pShell == 0)
	{
		m_pShell = new Console;
	}

	return m_pShell;
}

//--------------------------------------------------------------------
//! Sets the title of the console window.
void Console::SetTitle(const char* sz, ...)
{
	if (m_bActive)
	{
		// get a pointer to the argument list
		va_list	args;

		// make the message
		char sztitle[512];
		va_start(args, sz);
		vsprintf(sztitle, sz, args);
		va_end(args);

#ifdef WIN32
		SetConsoleTitleA(sztitle);
#endif
#ifdef LINUX
		printf("%c]0;%s%c", '\033', sztitle, '\007');
#endif
	}
}

//--------------------------------------------------------------------
//! get a command from the user

void Console::GetCommand(int& nargs, char **argv)
{
	static char szcmd[512] = {0};

	// print the command prompt
	printf("\n>");

	// you must flush the input buffer before using gets
	fflush(stdin);

	// get the command
	fgets(szcmd, 255, stdin);

	// fgets does not remove '\n' so we'll do it ourselves
	char* ch = strrchr(szcmd, '\n');
	if (ch) *ch = 0;

	// parse the arguments
	nargs = 0;
	int n = 0;
	int l = strlen(szcmd);
	ch = szcmd;
	for (int i=0; i<=l; ++i, ++ch)
	{
		if (!isspace(*ch) && (*ch != 0))
		{
			if (n==0) argv[nargs++] = ch;
			++n;
		}
		else 
		{
			if (n!=0)
			{
				argv[nargs-1][n] = 0;
				n = 0;
			}
		}
	}
}
