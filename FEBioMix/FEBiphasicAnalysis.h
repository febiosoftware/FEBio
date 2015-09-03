#pragma once
#include "FECore/FEAnalysis.h"
using namespace FECore;

//-----------------------------------------------------------------------------
//! Analysis class for biphasic problems
class FEBiphasicAnalysis : public FEAnalysis
{
public:
	FEBiphasicAnalysis(FEModel* pfem);
};
