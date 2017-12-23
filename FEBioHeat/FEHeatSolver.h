#pragma once

#include "FECore/FELinearSolver.h"
#include "FEHeatSolidDomain.h"

//-----------------------------------------------------------------------------
//! The FEHeatSolver solves heat-conduction problems
//! 

class FEHeatSolver : public FELinearSolver
{
public:
	//! constructor
	FEHeatSolver(FEModel* pfem);

	//! destructor
	~FEHeatSolver();

	//! one-time initialization
	bool Init() override;

protected: // from FELinearSolver

	//! calculate the residual
	void ForceVector(FEGlobalVector& R) override;

	//! calculate the stiffness matrix
	bool StiffnessMatrix(FELinearSystem& LS) override; 

protected:	// RHS helper functions
	//! Nodal fluxes
	void NodalFluxes(FEGlobalVector& R);

	//! Surface fluxes
	void SurfaceFluxes(FEGlobalVector& R);

	//! Heat Sources
	void HeatSources(FEGlobalVector& R);

protected:
	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};
