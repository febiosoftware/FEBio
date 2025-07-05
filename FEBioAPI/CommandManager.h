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



#pragma once

#include "Command.h"
#include <list>

class CommandManager
{
public:
	static CommandManager* GetInstance()
	{ 
		static bool bfirst = true;
		if (bfirst) { m_pMngr = new CommandManager; bfirst = false; }
		return m_pMngr; 
	}

public:
	void AddCommand(Command* pcmd) { m_Cmd.push_back(pcmd); }
	int Size() { return (int)m_Cmd.size(); }

	Command* Find(const char* szcmd)
	{
		int N = (int)m_Cmd.size();
		if (N == 0) return 0;

		std::list<Command*>::iterator ic = m_Cmd.begin();
		for (int i=0; i<N; ++i, ++ic)
		{
			if (strcmp(szcmd, (*ic)->GetName()) == 0) return (*ic);
		}

		return 0;
	}

	typedef std::list<Command*>::iterator CmdIterator;

	CmdIterator First() { return m_Cmd.begin(); }

protected:
	std::list<Command*>	m_Cmd;

protected:
	static CommandManager* m_pMngr;
};
