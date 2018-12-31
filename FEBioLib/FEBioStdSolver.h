#pragma once
#include <FECore/FECoreTask.h>
#include <FECore/FECoreKernel.h>

//-----------------------------------------------------------------------------
// This is the most commenly used task which will run a user-specified input 
// file. The results are stored in the logfile and the plotfile.
class FEBioStdSolver : public FECoreTask
{
public:
	FEBioStdSolver(FEModel* pfem);

	//! initialization
	bool Init(const char* szfile);

	//! Run the FE model
	bool Run();
};

//-----------------------------------------------------------------------------
class FEBioRestart : public FECoreTask
{
public:
	FEBioRestart(FEModel* pfem) : FECoreTask(pfem){}

	//! initialization
	bool Init(const char* szfile);

	//! Run the FE model
	virtual bool Run();
};
