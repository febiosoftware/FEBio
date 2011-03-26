#pragma once
#include "FESurfaceLoad.h"

//-----------------------------------------------------------------------------
//! Surface that sustains a heat flux boundary condition
//!
class FEHeatFlux : public FESurfaceLoad
{
public:
	struct LOAD
	{
		double	s[4];		// nodal scale factors
		int		lc;			// load curve
		LOAD()  { s[0] = s[1] = s[2] = s[3] = 1.0; lc = -1; }
	};

public:
	FEHeatFlux(FEMesh* pm) : FESurfaceLoad(pm){}

	//! allocate storage
	void create(int n)
	{
		m_surf.create(n);
		m_FC.resize(n);
	}

	//! clone
/*	FEDomain* Clone()
	{
		FEHeatFlux* ps = new FEHeatFlux(m_pMesh);
		ps->m_FC = m_FC;
		return ps;
	}
*/

	//! get a heat flux load BC
	LOAD& HeatFlux(int n) { return m_FC[n]; }

	//! stiffness matrix
	void StiffnessMatrix(FESolver* psolver) {}
	
	//! residual
	void Residual(FESolver* psolver, vector<double>& R);

protected:
	vector<LOAD>	m_FC;
};
