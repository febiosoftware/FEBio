#pragma once
#include <FECore/FEBioTask.h>

//-----------------------------------------------------------------------------
// This class defines the FEBioOpt task. It allows FEBio to call the optimization
// module to solve a parameter optimization problem.
class FEBioOpt : public FEBioTask
{
public:
	//! class constructor
	FEBioOpt(FEModel* pfem);

	//! Run the optimization module
	bool Run(const char* szfile);
};
