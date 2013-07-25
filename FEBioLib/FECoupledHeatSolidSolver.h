#pragma once
#include <FECore/FESolver.h>
#include <FEBioHeat/FEHeatSolver.h>
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

	//! Clean
	void Clean();

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

private: // not used
	virtual void AssembleResidual(vector<int>& en, vector<int>& elm, vector<double>& fe, vector<double>& R) { assert(false); }
	virtual void AssembleStiffness(vector<int>& en, vector<int>& elm, matrix& ke) { assert(false); }

protected:
	FEHeatSolver		m_Heat;		//!< heat solver
	FELinearSolidSolver	m_Solid;	//!< linear solid solver
};
