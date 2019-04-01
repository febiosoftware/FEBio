#include "stdafx.h"
#include "FEReactionRateConst.h"

// Material parameters for the FEReactionRateConst material
BEGIN_FECORE_CLASS(FEReactionRateConst, FEReactionRate)
	ADD_PARAMETER(m_k, FE_RANGE_GREATER_OR_EQUAL(0.0), "k");
END_FECORE_CLASS();
