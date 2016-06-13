#pragma once

#include "Command.h"
#include <list>
using namespace std;

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

		list<Command*>::iterator ic = m_Cmd.begin();
		for (int i=0; i<N; ++i, ++ic)
		{
			if (strcmp(szcmd, (*ic)->GetName()) == 0) return (*ic);
		}

		return 0;
	}

	typedef list<Command*>::iterator CmdIterator;

	CmdIterator First() { return m_Cmd.begin(); }

protected:
	list<Command*>	m_Cmd;

protected:
	static CommandManager* m_pMngr;
};
