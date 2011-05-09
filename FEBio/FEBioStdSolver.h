#pragma once
#include "FECore/FEBioTask.h"

//-----------------------------------------------------------------------------
// This is the most commenly used task which will run a user-specified input 
// file. The results are stored in the logfile and the plotfile.
class FEBioStdSolver : public FEBioTask
{
public:
	//! Run the FE model
	virtual bool Run(FEModel& fem, const char* szfile);
};

//-----------------------------------------------------------------------------
class FEBioStdSolverFactory : public FEBioTaskFactory
{
public:
	FEBioTask* CreateTask() { return new FEBioStdSolver(); }
};


//-----------------------------------------------------------------------------
class FEBioOptimize : public FEBioTask
{
public:
	//! Run the FE model
	virtual bool Run(FEModel& fem, const char* szfile);
};

//-----------------------------------------------------------------------------
class FEBioOptimizeFactory : public FEBioTaskFactory
{
public:
	FEBioTask* CreateTask() { return new FEBioOptimize(); }
};

//-----------------------------------------------------------------------------
class FEBioDiagnostic : public FEBioTask
{
public:
	//! Run the FE model
	virtual bool Run(FEModel& fem, const char* szfile);
};

//-----------------------------------------------------------------------------
class FEBioDiagnosticFactory : public FEBioTaskFactory
{
public:
	FEBioTask* CreateTask() { return new FEBioDiagnostic(); }
};


//-----------------------------------------------------------------------------
class FEBioRestart : public FEBioTask
{
public:
	//! Run the FE model
	virtual bool Run(FEModel& fem, const char* szfile);
};

//-----------------------------------------------------------------------------
class FEBioRestartFactory : public FEBioTaskFactory
{
public:
	FEBioTask* CreateTask() { return new FEBioRestart(); }
};
