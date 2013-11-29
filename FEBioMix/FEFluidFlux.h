#pragma once
#include "FECore/FESurfaceLoad.h"

//-----------------------------------------------------------------------------
//! This boundary condition sustains a fluid flux on a surface
//!
class FEFluidFlux : public FESurfaceLoad
{
public:
	struct LOAD
	{
		double	s[8];		// nodal scale factors
		int		lc;			// load curve
		int		bc;			// degree of freedom

		LOAD() { s[0] = s[1] = s[2] = s[3] = s[4] = s[5] = s[6] = s[7] = 1.0; bc = 0; }
	};

public:
	//! constructor
	FEFluidFlux(FEModel* pfem) : FESurfaceLoad(pfem) { m_blinear = false; m_bmixture = false; }

	void SetLinear(bool blinear) { m_blinear = blinear; }

	void SetMixture(bool bmix) { m_bmixture = bmix; }

	//! allocate storage
	void create(int n) { m_PC.resize(n); }

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
	void Residual(FEGlobalVector& R);

	//! serialize data
	void Serialize(DumpFile& ar);

protected:
	//! calculate stiffness for an element
	void FluxStiffness(FESurfaceElement& el, matrix& ke, vector<double>& vn, double dt, bool mixture);

	//! calculate stiffness for an element, for steady-state analysis
	void FluxStiffnessSS(FESurfaceElement& el, matrix& ke, vector<double>& vn, double dt, bool mixture);

	//! Calculates volumetric flow rate due to flux
	bool FlowRate(FESurfaceElement& el, vector<double>& fe, vector<double>& vn, double dt, bool mixture);

	//! Calculates volumetric flow rate due to flux, for steady-state analysis
	bool FlowRateSS(FESurfaceElement& el, vector<double>& fe, vector<double>& vn, double dt, bool mixture);

	//! Calculates the linear volumetric flow rate due to flux (ie. non-follower)
	bool LinearFlowRate(FESurfaceElement& el, vector<double>& fe, vector<double>& vn, double dt, bool mixture);

	//! Calculates the linear volumetric flow rate due to flux (ie. non-follower), for steady-state analysis
	bool LinearFlowRateSS(FESurfaceElement& el, vector<double>& fe, vector<double>& vn, double dt, bool mixture);

protected:
	bool	m_bmixture;		//!< mixture velocity or relative fluid flux
	bool	m_blinear;		//!< type (linear or nonlinear)

	// Fluid flux boundary data
	vector<LOAD>	m_PC;		//!< fluid flux boundary cards
};
