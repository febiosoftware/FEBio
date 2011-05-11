#pragma once

//-----------------------------------------------------------------------------
class FEModel;
class FEBioTask;

//-----------------------------------------------------------------------------
// The FEBioTask class is the base class for all tasks that FEBio can execute.
// A task is simply the highest level module which defines what FEBio will do.
// For example, the FEBioStdSolver will run a user-specified input file; the
// FEBioNAGOptimize does a material parameter optimization using the NAG library;
// the FEBioMaterialDiagnostic runs a material diagnostic on a set of predefined
// scenarios, etc. 
class FEBioTask
{
public:
	FEBioTask(void){}
	virtual ~FEBioTask(void){}

	//! Run the task.
	virtual bool Run (FEModel& fem, const char* szfile) = 0;
};
