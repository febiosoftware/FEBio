#pragma once
#include "FESurface.h"
#include "FEBoundaryCondition.h"

//-----------------------------------------------------------------------------
//! This class describes a fluid flux on a surface element
//! bc = 0 for solvent flux

class FEFluidFlux : public FEBoundaryCondition
{
public:
	FEFluidFlux() { s[0] = s[1] = s[2] = s[3] = 1.0; bc = 0; blinear = false; mixture = false; }

public:
	double	s[4];		// nodal scale factors
	int		face;		// face number
	int		lc;			// load curve
	int		bc;			// degree of freedom
	bool	blinear;	// linear or not (true is non-follower, false is follower)
	bool	mixture;	// mixture velocity or relative fluid flux
};

//-----------------------------------------------------------------------------
//! The flux surface is a surface domain that sustains a fluid flux boundary
//! condition
//!
class FEFluxSurface : public FESurface
{
public:
	//! constructor
	FEFluxSurface(FEMesh* pm) : FESurface(pm) {}

	//! allocate storage
	void create(int n)
	{
		FESurface::create(n);
		m_PC.resize(n);
	}

	//! clone
	FEDomain* Clone()
	{
		FEFluxSurface* ps = new FEFluxSurface(m_pMesh);
		ps->m_PC = m_PC;
		return ps;
	}

	//! get a flux BC
	FEFluidFlux& FluidFlux(int n) { return m_PC[n]; }

	//! calculate flux stiffness
	void StiffnessMatrix(FESolidSolver* psolver);

	//! calculate residual
	void Residual(FESolidSolver* psolver, vector<double>& R);

	//! serialize data
	void Serialize(FEM& fem, DumpFile& ar);

protected:
	//! calculate stiffness for an element
	void FluxStiffness(FESurfaceElement& el, matrix& ke, vector<double>& vn, double dt, bool mixture);

	//! Calculates volumetric flow rate due to flux
	bool FlowRate(FESurfaceElement& el, vector<double>& fe, vector<double>& vn, double dt, bool mixture);

	//! Calculates the linear volumetric flow rate due to flux (ie. non-follower)
	bool LinearFlowRate(FESurfaceElement& el, vector<double>& fe, vector<double>& vn, double dt, bool mixture);

protected:
	// Fluid flux boundary data
	vector<FEFluidFlux>	m_PC;		//!< fluid flux boundary cards
};
