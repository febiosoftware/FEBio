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
	while (m_Task.empty() == false) RemoveTask(m_Task[0]);
}

//-----------------------------------------------------------------------------
CTask* TMSession::GetTask(int i)
{
	if ((i>=0) && (i<(int)m_Task.size())) return m_Task[i];
	else return 0;
}

//-----------------------------------------------------------------------------
CTask* TMSession::GetVisibleTask(int n)
{
	if (n<0) return 0;
	int c = 0;
	for (int i=0; i<(int)m_Task.size(); ++i)
	{
		CTask* pt = m_Task[i];
		if (pt->IsVisible())
		{
			if (c == n) return pt; else c++;
		}
	}
	return 0;
}

//-----------------------------------------------------------------------------
CTask* TMSession::AddTask(const char* szfile)
{
	// see if this task already exists
	for (int i=0; i<(int)m_Task.size(); ++i)
	{
		CTask* pt = m_Task[i];
		if (strcmp(pt->GetFileName(), szfile) == 0) return 0;
	}

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
int TMSession::VisibleTasks()
{
	int nv = 0;
	for (int i=0; i<Tasks(); ++i) if (m_Task[i]->IsVisible()) nv++;
	return nv;
}

//-----------------------------------------------------------------------------
void TMSession::RemoveTask(CTask* pt)
{
	assert(pt);
	for (int i=0; i<(int)m_Task.size(); ++i)
	{
		CTask* pti = GetTask(i);
		if (pt == pti)
		{
			vector<CTask*>::iterator it = m_Task.begin() + i;
			m_Task.erase(it);
			delete pti;
			break;
		}
	}
}
