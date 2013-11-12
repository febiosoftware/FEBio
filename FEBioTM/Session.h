#pragma once
#include "Task.h"
#include <vector>
using namespace std;

//-----------------------------------------------------------------------------
// This class defines an FEBioTM session
class TMSession
{
public:
	TMSession();
	~TMSession();

	// clear the session
	void Clear();

	// add a taks to the document
	CTask* AddTask(const char* szfile);

	// remove a task from the queue
	void RemoveTask(int n);

	// get the number of tasks
	int Tasks() { return (int) m_Task.size(); }

	// get a task
	CTask* GetTask(int i);

private:
	vector<CTask*>	m_Task;
};
