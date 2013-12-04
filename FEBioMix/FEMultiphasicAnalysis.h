#pragma once
#include "FECore/FEAnalysis.h"
using namespace FECore;

//-----------------------------------------------------------------------------
//! Analysis class for multiphasic problems
class FEMultiphasicAnalysis : public FEAnalysis
{
public:
	FEMultiphasicAnalysis(FEModel* pfem) : FEAnalysis(pfem, FE_MULTIPHASIC) {}

	bool Init();

protected:
	void InitNodes();
};
