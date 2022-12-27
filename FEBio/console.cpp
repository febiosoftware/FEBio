/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
the City of New York, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/



#include "stdafx.h"
#include "console.h"
#include <stdarg.h>
#include <string.h>
#include <ctype.h>

#ifdef WIN32
	#include <shobjidl.h>
	#include <windows.h>
	#include <conio.h>
#else
	//These are for the wait command
	#include <termios.h>
	#include <unistd.h>
#endif

#ifdef WIN32
ITaskbarList3*	taskBar = NULL;
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

Console::~Console()
{
	CleanUp();
}

void Console::CleanUp()
{
#ifdef WIN32
	if (taskBar != NULL)
	{
		taskBar->Release();
		CoUninitialize();
		taskBar = NULL;
	}
#endif
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
//! gets a line from the user
void Console::Wait()
{
	// notify the user.
	fprintf(stderr, "Press any key to continue...\n");

	// flush the standard stream
	fflush(stdin);

	// wait until user inputs a character with no return symbol.
#ifdef WIN32
	_getch();
#else
	// change mode of termios so echo is off.
	struct termios oldt,
					newt;
	tcgetattr( STDIN_FILENO, &oldt );
	newt = oldt;
	newt.c_lflag &= ~( ICANON | ECHO );
	tcsetattr( STDIN_FILENO, TCSANOW, &newt );
	getchar();
	// reset stdin.
	tcsetattr( STDIN_FILENO, TCSANOW, &oldt );
#endif

}

//--------------------------------------------------------------------
//! get a command from the user

void Console::GetCommand(int& nargs, char **argv)
{
	static char szcmd[512] = { 0 };
	szcmd[0] = 0;

	// print the command prompt
	Write("\nfebio>", 0x0E);

	// you must flush the input buffer before using gets
	fflush(stdin);

	// get the command
	fgets(szcmd, 255, stdin);

	// fgets does not remove '\n' so we'll do it ourselves
	char* ch = strrchr(szcmd, '\n');
	if (ch) *ch = 0;

	// check for a percentage sign
	if (szcmd[0] == '%')
	{
		int n = atoi(szcmd + 1) - 1;
		if ((n >= 0) && (n < m_history.size()))
		{
			strcpy(szcmd, m_history[n].c_str());
		}
		else { nargs = 0; return; }
	}

	// store a copy of the input to the history
	// (unless it's the history command)
	if (strcmp(szcmd, "hist") != 0)
		m_history.push_back(szcmd);

	// parse the arguments
	nargs = 0;
	int n = 0;
	int b = 0;
	ch = szcmd;
	while (*ch)
	{
		switch (*ch)
		{
		case ' ':
			if ((b == 0) && (n != 0))
			{
				*ch = 0;
				n = 0;
			}
			break;
		case '"':
			if ((b == 0) && (n==0)) b = 1;
			else 
			{
				b = 0;
				*ch = 0;
				n = 0;
			}
			break;
		default:
			if (n == 0) argv[nargs++] = ch;
			n++;
		}
		ch++;
	}
}

//--------------------------------------------------------------------
const std::vector<std::string>& Console::GetHistory() const
{
	return m_history;
}

//--------------------------------------------------------------------
//! this function draws an image to the console

void Console::Draw(unsigned char *img, int nx, int ny)
{
#ifdef WIN32
	int i, j, n = 0;
	WORD att;
	printf("\n");
	HANDLE hout = GetStdHandle(STD_OUTPUT_HANDLE);
	DWORD col[] = {0x00, 0x04, 0x02, 0x01, 0x0C, 0x0A, 0x09, 0x08, 0x07};
	for (j=0; j<ny; ++j)
	{
		for (i=0; i<nx; ++i)
		{
			att = (WORD) ((col[img[n++]] << 4)%0xFF);
			SetConsoleTextAttribute(hout, att);
			printf(" ");
		}
		printf("\n");
	}
	SetConsoleTextAttribute(hout, 0x0F);
#endif
}

//--------------------------------------------------------------------
void Console::Write(const char *sz, unsigned short att)
{
#ifdef WIN32
	printf("\n");
	HANDLE hout = GetStdHandle(STD_OUTPUT_HANDLE);
	SetConsoleTextAttribute(hout, (WORD) att);
	printf("%s", sz);
	SetConsoleTextAttribute(hout, 0x0F);
#else
	printf("%s", sz);
#endif
}

//--------------------------------------------------------------------
void Console::SetProgress(double pct)
{
#ifdef WIN32
	// Get the console window's handle
	HWND hwnd = GetConsoleWindow();
	if (hwnd == NULL) return;

	// initialize task bar
	if (taskBar == NULL)
	{
		CoInitialize(NULL);
		CoCreateInstance(CLSID_TaskbarList, NULL, CLSCTX_INPROC_SERVER, IID_ITaskbarList3, (void **)&taskBar);
	}

	if (taskBar)
	{
		if ((pct <= 0.0) || (pct >= 100.0))
		{
			taskBar->SetProgressState(hwnd, TBPF_NOPROGRESS);
		}
		else
		{
			taskBar->SetProgressState(hwnd, TBPF_NORMAL);
			taskBar->SetProgressValue(hwnd, (ULONGLONG)pct, 100);
		}
	}
	else CoUninitialize();
#endif
}

void ConsoleStream::print(const char* sz)
{ 
	fprintf(stdout, "%s", sz);
	fflush(stdout);
}
