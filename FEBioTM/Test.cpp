#include "stdafx.h"
#include "Test.h"
#include <assert.h>

//-----------------------------------------------------------------------------
CTest::CTest(TMSession& session) : m_session(session)
{
}

//-----------------------------------------------------------------------------
CTest::~CTest()
{
}

//-----------------------------------------------------------------------------
bool CTest::LoadData(const char* szfile)
{
	FILE* fp = fopen(szfile, "rt");
	if (fp == 0) return false;

	m_data.clear();
	char szline[1025] = {0};
	while (fgets(szline, 1024, fp))
	{
		STATS s;
		int nread = sscanf(szline, "%s%d%d%d%d%d%d", s.szfile, &s.nreturn, &s.ntime, &s.niters, &s.nrhs, &s.nreform, &s.nsec);
		if (nread == 7) m_data.push_back(s);
	}

	int N = m_data.size();
	m_result.assign(N, 0);

	fclose(fp);
	return true;
}

//-----------------------------------------------------------------------------
bool CTest::SaveData(const char* szfile)
{
	FILE* fp = fopen(szfile, "wt");
	if (fp == 0) return false;

	int N = m_session.Tasks();
	for (int i=0; i<N; ++i)
	{
		CTask& task = *m_session.GetTask(i);
		CTask::STATS& s = task.m_stats;
		fprintf(fp, "%s %d %d %d %d %d %d\n", task.GetFileTitle(), s.nreturn, s.ntime, s.niters, s.nrhs, s.nreform, s.nsec);
	}

	fclose(fp);
	return true;
}

//-----------------------------------------------------------------------------
void CTest::Run()
{
	if (m_data.empty())
	{
		// collect data
		int N = m_session.Tasks();
		m_data.resize(N);
		m_result.assign(N, 2);
		for (int i=0; i<N; ++i)
		{
			CTask& task = *m_session.GetTask(i);
			STATS& s = m_data[i];
			strcpy(s.szfile, task.GetFileTitle());
			s.nreturn = task.m_stats.nreturn;
			s.ntime   = task.m_stats.ntime;
			s.niters  = task.m_stats.niters;
			s.nrhs    = task.m_stats.nrhs;
			s.nreform = task.m_stats.nreform;
			s.nsec    = task.m_stats.nsec;
		}
	}
	else
	{
		// compare results
		int N = m_session.Tasks();
		for (int i=0; i<N; ++i)
		{
			CTask& task = *m_session.GetTask(i);

			int n = FindStatsIndex(task.GetFileTitle());
			if (n == -1)
			{
				STATS s;
				strcpy(s.szfile, task.GetFileTitle());
				s.nreturn = task.m_stats.nreturn;
				s.ntime   = task.m_stats.ntime;
				s.niters  = task.m_stats.niters;
				s.nrhs    = task.m_stats.nrhs;
				s.nreform = task.m_stats.nreform;
				s.nsec    = task.m_stats.nsec;
				m_data.push_back(s);
				m_result.push_back(2);
			}
			else
			{
				m_result[n] = 1;
				if (m_data[n].nreturn != 1) m_result[n] = 0;
				if (m_data[n].nreturn != task.m_stats.nreturn) m_result[n] = 0;
				if (m_data[n].ntime   != task.m_stats.ntime  ) m_result[n] = 0;
				if (m_data[n].niters  != task.m_stats.niters ) m_result[n] = 0;
				if (m_data[n].nrhs    != task.m_stats.nrhs   ) m_result[n] = 0;
				if (m_data[n].nreform != task.m_stats.nreform) m_result[n] = 0;
			}
		}
	}
}

//-----------------------------------------------------------------------------
CTest::STATS CTest::GetStats(CTask& task)
{
	int n = FindStatsIndex(task.GetFileTitle());
	assert(n >= 0);
	return GetStats(n);
}

//-----------------------------------------------------------------------------
int CTest::FindStatsIndex(const char* szfile)
{
	int N = m_data.size();
	for (int i=0; i<N; ++i)
	{
		if (strcmp(m_data[i].szfile, szfile) == 0) return i;
	}
	return -1;
}
