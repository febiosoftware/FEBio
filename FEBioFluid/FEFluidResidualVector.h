#pragma once
#include "FECore/FEGlobalVector.h"
#include <vector>

//-----------------------------------------------------------------------------
class FEModel;

//-----------------------------------------------------------------------------
//! The FEFluidResidualVector implements a global vector that stores the residual.

class FEFluidResidualVector : public FEGlobalVector
{
public:
	//! constructor
	FEFluidResidualVector(FEModel& fem, std::vector<double>& R, std::vector<double>& Fr);

	//! destructor
	~FEFluidResidualVector();

	//! Assemble the element vector into this global vector
	void Assemble(vector<int>& en, vector<int>& elm, vector<double>& fe);
};
