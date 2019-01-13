#pragma once
#include "FECore/FEGlobalVector.h"
#include <vector>

//-----------------------------------------------------------------------------
class FEModel;

//-----------------------------------------------------------------------------
//! The FEResidualVector implements a global vector that stores the residual.

class FECORE_API FEResidualVector : public FEGlobalVector
{
public:
	//! constructor
	FEResidualVector(FEModel& fem, std::vector<double>& R, std::vector<double>& Fr);

	//! destructor
	~FEResidualVector();

	//! Assemble the element vector into this global vector
	void Assemble(vector<int>& en, vector<int>& elm, vector<double>& fe, bool bdom = false);
};
