#pragma once

#include <FEBioLib\FESolver.h>
#include "FEHeatSolidDomain.h"

//-----------------------------------------------------------------------------
//! The FEHeatSolver solves heat-conduction problems
//! 

class FEHeatSolver : public FESolver
{
public:
	//! constructor
	FEHeatSolver(FEModel& fem);

	//! one-time initialization
	bool Init();

	//! solve a timestep
	bool SolveStep(double time);

	//! serialize data
	void Serialize(DumpFile& ar);

	//! Initialize linear equation system
	bool InitEquations();

protected:
	//! calculate the residual
	void Residual();

	//! calculate the stiffness matrix
	bool StiffnessMatrix(); 

	//! form the stiffness matrix
	bool ReformStiffness();

	//! update solution
	void Update(vector<double>& u);

	//! Prep the step
	void PrepStep();

protected:	// Residual functions
	//! Nodal fluxes
	void NodalFluxes(FEGlobalVector& R);

	//! Surface fluxes
	void SurfaceFluxes(FEGlobalVector& R);

	//! Heat Sources
	void HeatSources(FEGlobalVector& R);

public:
	//! assemble element stiffness matrix
	void AssembleStiffness(vector<int>& en, vector<int>& lm, matrix& ke);

	//! assemble an element stiffness into the residual
	void AssembleResidual(vector<int>& lm, matrix& kc);

private: // TODO: use this function instead
//	virtual void AssembleResidual(vector<int>& en, vector<int>& elm, vector<double>& fe, vector<double>& R) { assert(false); }

public:
	vector<double>	m_T;	//!< temperature array
	vector<double>	m_Tp;	//!< previous temperatures
	vector<double>	m_R;	//!< residual
	vector<double>	m_u;	//!< prescribed temperatures

protected:
	bool	m_brhs;	//!< flag used to indicate if element stiffness must be assembled into RHS

	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};
