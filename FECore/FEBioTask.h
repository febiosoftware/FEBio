#pragma once
#include "FECoreBase.h"

//-----------------------------------------------------------------------------
class FEModel;

//-----------------------------------------------------------------------------
// The FEBioTask class is the base class for all tasks that FEBio can execute.
// A task is simply the highest level module which defines what FEBio will do.
// For example, the FEBioStdSolver will run a user-specified input file; the
// FEBioNAGOptimize does a material parameter optimization using the NAG library;
// the FEBioMaterialDiagnostic runs a material diagnostic on a set of predefined
// scenarios, etc. 
class FEBioTask : public FECoreBase
{
public:
	FEBioTask(FEModel* pfem) : FECoreBase(FETASK_ID) { m_pfem = pfem; }
	virtual ~FEBioTask(void){}

	//! Run the task.
	virtual bool Run (const char* szfile) = 0;

protected:
	FEModel*	m_pfem;
};
