#pragma once
#include "FEBiphasicSolver.h"

//-----------------------------------------------------------------------------
// This class adds additional functionality to the FESolidSolver to solve
// solute problems. 
class FEBiphasicSoluteSolver : public FEBiphasicSolver
{
public:
	//! con/descructor
	FEBiphasicSoluteSolver(FEModel& fem);
	virtual ~FEBiphasicSoluteSolver(){}

	//! Initialize data structures
	bool Init();

	//! prepares the data for the first QN iteration
	virtual void PrepStep(double time);

	//! Performs a Newton-Raphson iteration
	bool Quasin(double time);

	//! serialize data to/from dump file
	void Serialize(DumpFile& ar);

public:
	//! Calculates residual (overridden from FEBiphasicSolver)
	virtual bool Residual(vector<double>& R);

protected:
	void GetConcentrationData(vector<double>& ci, vector<double>& ui);

public:
	double	m_Ctol;			//!< concentration tolerance
};
