#pragma once
#include "FECore/FESolver.h"
#include "FEHeatSolver.h"
#include "FEBioMech/FELinearSolidSolver.h"

//-----------------------------------------------------------------------------
//! This class implements a coupled thermo-elastic solver
//!
class FECoupledHeatSolidSolver : public FESolver
{
public:
	//! constructor
	FECoupledHeatSolidSolver(FEModel* pfem);

	//! destructor
	~FECoupledHeatSolidSolver(){}

	//! Initializiation
	bool Init() override;

	//! Clean
	void Clean() override;

	//! Solve a step
	bool SolveStep() override;

	//! Update solution
	void Update(vector<double>& u) override;

	//! data serialization
	void Serialize(DumpStream& ar) override;

	//! Initialize equations
	bool InitEquations() override;

protected:
	//! calculate "initial" stresses base on temperatures
	void CalculateInitialStresses();

private: // not used
	virtual void AssembleResidual(vector<int>& en, vector<int>& elm, vector<double>& fe, vector<double>& R) { assert(false); }
	virtual void AssembleStiffness(vector<int>& en, vector<int>& elm, matrix& ke) override { assert(false); }

protected:
	FEHeatSolver		m_Heat;		//!< heat solver
	FELinearSolidSolver	m_Solid;	//!< linear solid solver
};
