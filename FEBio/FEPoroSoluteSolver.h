#pragma once
#include "FEPoroSolidSolver.h"

//-----------------------------------------------------------------------------
// This class adds additional functionality to the FESolidSolver to solve
// solute problems. 
class FEPoroSoluteSolver : public FEPoroSolidSolver
{
public:
	//! con/descructor
	FEPoroSoluteSolver(FEM& fem);
	virtual ~FEPoroSoluteSolver(){}

	//! Initialize data structures
	bool Init();

	//! prepares the data for the first QN iteration
	virtual void PrepStep(double time);

	//! Performs a Newton-Raphson iteration
	bool Quasin(double time);

	//! serialize data to/from dump file
	void Serialize(DumpFile& ar);

protected:
	void GetConcentrationData(vector<double>& ci, vector<double>& ui);

public:
	double	m_Ctol;			//!< concentration tolerance
};
