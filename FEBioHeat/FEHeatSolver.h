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
	bool Init();

	//! serialize data
	void Serialize(DumpFile& ar);

protected: // from FELinearSolver

	//! calculate the residual
	void RHSVector(FEGlobalVector& R);

	//! calculate the stiffness matrix
	bool StiffnessMatrix(); 

	//! update solution
	void Update(vector<double>& u);

protected:	// RHS helper functions
	//! Nodal fluxes
	void NodalFluxes(FEGlobalVector& R);

	//! Surface fluxes
	void SurfaceFluxes(FEGlobalVector& R);

	//! Heat Sources
	void HeatSources(FEGlobalVector& R);

public:
	//! assemble element stiffness matrix
	void AssembleStiffness(vector<int>& en, vector<int>& lm, matrix& ke);

public:
	vector<double>	m_Tp;	//!< previous temperatures

protected:
	bool	m_brhs;	//!< flag used to indicate if element stiffness must be assembled into RHS

	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};
