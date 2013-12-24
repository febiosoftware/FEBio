#pragma once
#include "FECore/FECoreTask.h"

//-----------------------------------------------------------------------------
// This class defines the FEBioOpt task. It allows FEBio to call the optimization
// module to solve a parameter optimization problem.
class FEBioOpt : public FECoreTask
{
public:
	//! class constructor
	FEBioOpt(FEModel* pfem);

	//! Run the optimization module
	bool Run(const char* szfile);
};
