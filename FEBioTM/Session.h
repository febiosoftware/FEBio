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
	void RemoveTask(CTask* pt);

	// get the number of tasks
	int Tasks() { return (int) m_Task.size(); }

	// return the number of visible tasks
	int VisibleTasks();

	// get a task
	CTask* GetTask(int i);

	// get the visible task
	CTask* GetVisibleTask(int i);

private:
	vector<CTask*>	m_Task;
};
