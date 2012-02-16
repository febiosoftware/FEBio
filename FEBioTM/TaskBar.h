#pragma once

#include <Flx_Group.h>

class CTask;

class CTaskBar : public Flx_Group
{
public:
	CTaskBar(int x, int y, int w, CTask* pt);
	~CTaskBar() { m_pTask = 0; }

protected:
	CTask*	m_pTask;
};
