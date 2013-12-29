#pragma once
#include "Session.h"
#include <vector>
using namespace std;

//-----------------------------------------------------------------------------
//! This class represents a test
class CTest
{
public:
	struct STATS
	{
		char	szfile[512];	// file name
		int		nreturn;	// return value
		int		ntime;		// number of time steps
		int		niters;		// number of iterations
		int		nrhs;		// number of RHS evaluations
		int		nreform;	// number of reformations
		int		nsec;		// number of seconds runtime
	};

public:
	CTest(TMSession& session);
	virtual ~CTest();

	void Run();

	bool SaveData(const char* szfile);
	bool LoadData(const char* szfile);

	STATS GetStats(int n) { return m_data[n]; }
	int GetResult(int n) { return m_result[n]; }

	STATS GetStats(CTask& task);

	int FindStatsIndex(const char* szfile);

private:
	TMSession&		m_session;
	vector<int>		m_result;
	vector<STATS>	m_data;
};
