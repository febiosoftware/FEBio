#include "stdafx.h"
#include "FEOptimizeData.h"
#include "FELMOptimizeMethod.h"
#include "FEOptimizeInput.h"
#include <FECore/FECoreKernel.h>

//-----------------------------------------------------------------------------
FEOptimizeData::FEOptimizeData(FEModel& fem) : m_fem(fem)
{
	m_pSolver = 0;
	m_pTask = 0;
}

//-----------------------------------------------------------------------------
FEOptimizeData::~FEOptimizeData(void)
{
	delete m_pSolver;
}

//-----------------------------------------------------------------------------
bool FEOptimizeData::Init()
{
	// allocate default optimization solver if none specified in input file
	if (m_pSolver == 0) m_pSolver = new FELMOptimizeMethod;

	// allocate default solver if none specified in input file
	if (m_pTask == 0) m_pTask = fecore_new<FECoreTask>(FETASK_ID, "solve", &m_fem);

	// do the initialization of the task
	if (m_pTask->Init(0) == false) return false;

	// find the variable
	OPT_OBJECTIVE& obj = GetObjective();
	obj.m_pd = m_fem.FindParameter(obj.m_szname);
	if (obj.m_pd == 0) return false;

	return true;
}

//-----------------------------------------------------------------------------
bool FEOptimizeData::Solve()
{
	// make sure we have a task that will solve the FE model
	if (m_pTask == 0) return false;

	// go for it!
	return m_pSolver->Solve(this);
}

//-----------------------------------------------------------------------------
//! Read the data from the input file
//!
bool FEOptimizeData::Input(const char *szfile)
{
	FEOptimizeInput in;
	if (in.Input(szfile, this) == false) return false;
	return true;
}

//-----------------------------------------------------------------------------
bool FEOptimizeData::RunTask()
{
	if (m_pTask == 0) return false;
	return m_pTask->Run();
}
