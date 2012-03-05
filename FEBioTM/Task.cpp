#include "stdafx.h"
#include "Task.h"
#include <assert.h>

//-----------------------------------------------------------------------------
// initialize static CTask variables
CTask* CTask::m_prun = 0;

//-----------------------------------------------------------------------------
void CTask::SetFileName(const char* szfile)
{
	m_szfile[0] = 0;
	int l = strlen(szfile)+1;
	assert((l>1) && (l<MAX_FILE));
	if ((l > 1) && (l<MAX_FILE)) strncpy(m_szfile, szfile, l);
}

//-----------------------------------------------------------------------------
const char* CTask::GetFileTitle()
{
	char* c1 = strrchr(m_szfile, '\\');
	char* c2 = strrchr(m_szfile, '/');
	if ((c1 == 0) && (c2 == 0)) return m_szfile;
	if (c1 == 0) return c2+1;
	if (c2 == 0) return c1+1;
	if (c1 < c2) return c2+1; else return c1+1;
}

//-----------------------------------------------------------------------------
void CTask::GetFilePath(char* szpath)
{
	char* c1 = strrchr(m_szfile, '\\');
	char* c2 = strrchr(m_szfile, '/');
	if ((c1 == 0) && (c2 == 0)) strcpy(szpath, m_szfile);
	if ((c1 == 0) || (c2 > c1)) strncpy(szpath, m_szfile, c2 - m_szfile);
	if ((c2 == 0) || (c1 > c2)) strncpy(szpath, m_szfile, c1 - m_szfile);
}
