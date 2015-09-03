#pragma once
#include "FECore/FEAnalysis.h"
using namespace FECore;

//-----------------------------------------------------------------------------
//! This class describes a heat-transfer analysis
class FEHeatTransferAnalysis : public FEAnalysis
{
public:
	//! constructor
	FEHeatTransferAnalysis(FEModel* pfem);
};
