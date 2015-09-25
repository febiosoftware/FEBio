#pragma once
#include "FESolidSolver.h"

//-----------------------------------------------------------------------------
//! This class implements a solver for solid mechanics problems that uses
//! the conjugate gradient method to solve the nonlinear finite element equations
class FECGSolidSolver : public FESolidSolver
{
public:
	//! constructor
	FECGSolidSolver(FEModel* pfem);

	//! initialization
	bool Init();

	//! Performs a CG step
	bool Quasin(double time);

	//! update nodal positions, velocities, accelerations, etc.
	void UpdateKinematics(vector<double>& ui);

private:
	//! modified linesearch for Hager-Zhang solver
	double LineSearchCG(double s);

};
