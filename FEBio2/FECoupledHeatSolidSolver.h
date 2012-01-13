#pragma once
#include "FESolver.h"
#include "FEHeatSolver.h"
#include "FELinearSolidSolver.h"

//-----------------------------------------------------------------------------
//! This class implements a coupled thermo-elastic solver
//!
class FECoupledHeatSolidSolver : public FESolver
{
public:
	//! constructor
	FECoupledHeatSolidSolver(FEModel& fem);

	//! destructor
	~FECoupledHeatSolidSolver(){}

	//! Initializiation
	bool Init();

	//! Solve a step
	bool SolveStep(double time);

	//! Update solution
	void Update(vector<double>& u);

	//! data serialization
	void Serialize(DumpFile& ar);

	//! Initialize equations
	bool InitEquations();

protected:
	//! calculate "initial" stresses base on temperatures
	void CalculateInitialStresses();

protected:
	FEHeatSolver		m_Heat;		//!< heat solver
	FELinearSolidSolver	m_Solid;	//!< linear solid solver
};
