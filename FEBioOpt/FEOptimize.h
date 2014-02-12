#pragma once
#include "FECore/FECoreTask.h"

//-----------------------------------------------------------------------------
// This class defines the FEOptimize task. It allows FEBio to call the optimization
// module to solve a parameter optimization problem.
class FEOptimize : public FECoreTask
{
public:
	//! class constructor
	FEOptimize(FEModel* pfem);

	//! Run the optimization module
	bool Run(const char* szfile);
};
