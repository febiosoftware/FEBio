// console.cpp: implementation of the Console class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "console.h"
#include <stdarg.h>

#ifdef WIN32
	#include <windows.h>
#endif

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

bool Console::m_bActive = true;

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
