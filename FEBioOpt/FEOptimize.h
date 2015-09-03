#pragma once
#include "FECore/FECoreTask.h"
#include "FEOptimizer.h"

//-----------------------------------------------------------------------------
// This class defines the FEOptimize task. It allows FEBio to call the optimization
// module to solve a parameter optimization problem.
class FEOptimize : public FECoreTask
{
public:
	//! class constructor
	FEOptimize(FEModel* pfem);

	//! initialization
	bool Init(const char* szfile);

	//! Run the optimization module
	bool Run();

private:
	FEOptimizeData	m_opt;
};
