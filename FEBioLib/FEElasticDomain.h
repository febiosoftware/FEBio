#pragma once
#include "FECore/FENLSolver.h"
#include "FECore/FEModel.h"
#include <vector>
using namespace std;

//-----------------------------------------------------------------------------
//! Abstract interface class for elastic domains
class FEElasticDomain
{
public:
	virtual ~FEElasticDomain(){}
	virtual void StiffnessMatrix(FENLSolver* psolver) = 0;
	virtual void Residual(FENLSolver* psolver, vector<double>& R) = 0;
	virtual void UpdateStresses(FEModel& fem) = 0;
};
