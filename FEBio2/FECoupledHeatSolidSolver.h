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
	FECoupledHeatSolidSolver(FEM& fem);

	//! destructor
	~FECoupledHeatSolidSolver(){}

	//! Initializiation
	bool Init();

	//! Solve a step
	bool SolveStep(double time);

	//! data serialization
	void Serialize(DumpFile& ar);

protected:
	FEHeatSolver		m_Heat;		//!< heat solver
	FELinearSolidSolver	m_Solid;	//!< linear solid solver
};
