#pragma once
#include "FECore/FEBioTask.h"
#include "FECore/febio.h"

//-----------------------------------------------------------------------------
// This is the most commenly used task which will run a user-specified input 
// file. The results are stored in the logfile and the plotfile.
class FEBioStdSolver : public FEBioTask
{
public:
	FEBioStdSolver(FEModel* pfem) : FEBioTask(pfem){}

	//! Run the FE model
	virtual bool Run(const char* szfile);
};

//-----------------------------------------------------------------------------
class FEBioDiagnostic : public FEBioTask
{
public:
	FEBioDiagnostic(FEModel* pfem) : FEBioTask(pfem){}
	//! Run the FE model
	virtual bool Run(const char* szfile);
};

//-----------------------------------------------------------------------------
class FEBioRestart : public FEBioTask
{
public:
	FEBioRestart(FEModel* pfem) : FEBioTask(pfem){}
	//! Run the FE model
	virtual bool Run(const char* szfile);
};
