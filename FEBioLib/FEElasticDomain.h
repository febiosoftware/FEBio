#pragma once
#include "FECore/FENLSolver.h"
#include "FECore/FEModel.h"
#include "FECore/FEBodyForce.h"
#include <vector>
using namespace std;

//-----------------------------------------------------------------------------
//! Abstract interface class for elastic domains
class FEElasticDomain
{
public:
	virtual ~FEElasticDomain(){}
	virtual void StiffnessMatrix(FENLSolver* psolver) = 0;
	virtual void InertialForces(FENLSolver* psolver, vector<double>& R, vector<double>& F) = 0;
	virtual void Residual(FENLSolver* psolver, vector<double>& R) = 0;
	virtual void UpdateStresses(FEModel& fem) = 0;

	// these functions are replacing Residual
	virtual void InternalForces(FENLSolver* psolver, vector<double>& R) = 0;
	virtual void BodyForce(FENLSolver* psolver, FEBodyForce& bf, vector<double>& R) = 0;
};
