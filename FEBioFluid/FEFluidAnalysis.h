#pragma once
#include "FECore/FEAnalysis.h"
using namespace FECore;

//-----------------------------------------------------------------------------
//! Analysis class for fluid problems
class FEFluidAnalysis : public FEAnalysis
{
public:
	FEFluidAnalysis(FEModel* pfem);
};
