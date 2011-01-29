#pragma once
#include "FESolidSolver.h"

//-----------------------------------------------------------------------------
// This class adds additional functionality to the FESolidSolver to solve
// poro/solute problems. 
class FEPoroSolidSolver : public FESolidSolver
{
public:
	//! constructor
	FEPoroSolidSolver(FEM& fem) : FESolidSolver(fem) {}

	//! Performs a Newton-Raphson iteration
	bool Quasin(double time);
};
