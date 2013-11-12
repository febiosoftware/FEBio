#include "stdafx.h"
#include "Session.h"
#include <assert.h>

//-----------------------------------------------------------------------------
TMSession::TMSession()
{
}

//-----------------------------------------------------------------------------
TMSession::~TMSession()
{
	Clear();
}

//-----------------------------------------------------------------------------
void TMSession::Clear()
{
	while (m_Task.empty() == false) RemoveTask(0);
}

//-----------------------------------------------------------------------------
CTask* TMSession::GetTask(int i)
{
	if ((i>=0) && (i<(int)m_Task.size())) return m_Task[i];
	else return 0;
}

//-----------------------------------------------------------------------------
CTask* TMSession::AddTask(const char* szfile)
{
	// create a new task
	CTask* pt = new CTask;
	pt->SetFileName(szfile);

	// create a text buffer
	Fl_Text_Buffer* pb = new Fl_Text_Buffer;
	int nret = pb->appendfile(szfile);
	pt->SetTextBuffer(pb);

	m_Task.push_back(pt);
	return pt;
}

//-----------------------------------------------------------------------------
void TMSession::RemoveTask(int n)
{
	CTask* pt = GetTask(n);
	assert(pt);
	if (pt == 0) return;
	vector<CTask*>::iterator it = m_Task.begin() + n;
	m_Task.erase(it);
	delete pt;
}
