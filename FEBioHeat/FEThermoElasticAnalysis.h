#pragma once
#include <FECore/FEAnalysis.h>
using namespace FECore;

//-----------------------------------------------------------------------------
class FEThermoElasticAnalysis : public FEAnalysis
{
public:
	//! constructor
	FEThermoElasticAnalysis(FEModel* pfem);

	//! Initialization
	bool Init();
};
