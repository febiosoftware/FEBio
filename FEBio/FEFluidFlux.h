#pragma once
#include "FESurfaceLoad.h"

//-----------------------------------------------------------------------------
//! This boundary condition sustains a fluid flux on a surface
//!
class FEFluidFlux : public FESurfaceLoad
{
public:
	struct LOAD
	{
		double	s[4];		// nodal scale factors
		int		lc;			// load curve
		int		bc;			// degree of freedom

		LOAD() { s[0] = s[1] = s[2] = s[3] = 1.0; bc = 0; }
	};

public:
	//! constructor
	FEFluidFlux(FEMesh* pm, bool blinear = false, bool bmixture = false) : FESurfaceLoad(pm) { m_blinear = blinear; m_bmixture = bmixture; }

	//! allocate storage
	void create(int n)
	{
		m_surf.create(n);
		m_PC.resize(n);
	}

	//! clone
/*	FEDomain* Clone()
	{
		FEFluidFlux* ps = new FEFluidFlux(m_pMesh);
		ps->m_PC = m_PC;
		return ps;
	}
*/

	//! get a flux BC
	LOAD& FluidFlux(int n) { return m_PC[n]; }

	//! calculate flux stiffness
	void StiffnessMatrix(FESolver* psolver);

	//! calculate residual
	void Residual(FESolver* psolver, vector<double>& R);

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
	bool	m_bmixture;		//!< mixture velocity or relative fluid flux
	bool	m_blinear;		//!< type (linear or nonlinear)

	// Fluid flux boundary data
	vector<LOAD>	m_PC;		//!< fluid flux boundary cards
};
