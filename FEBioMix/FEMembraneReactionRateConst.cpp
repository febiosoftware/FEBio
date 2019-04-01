#include "stdafx.h"
#include "FEMembraneReactionRateConst.h"

// Material parameters for the FEMembraneReactionRateConst material
BEGIN_FECORE_CLASS(FEMembraneReactionRateConst, FEMembraneReactionRate)
	ADD_PARAMETER(m_k, FE_RANGE_GREATER_OR_EQUAL(0.0), "k");
END_FECORE_CLASS();
