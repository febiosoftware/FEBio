#pragma once
#include "FESurface.h"
#include "FEBoundaryCondition.h"

//-----------------------------------------------------------------------------
//! This class describes a heat flux card
//!
class FEHeatFlux : public FEBoundaryCondition
{
public:
	FEHeatFlux()  { s[0] = s[1] = s[2] = s[3] = 1.0; lc = -1; face = -1; }

public:
	double	s[4];		// nodal scale factors
	int		face;		// face number
	int		lc;			// load curve
};

//-----------------------------------------------------------------------------
//! Surface that sustains a heat flux boundary condition
//!
class FEHeatFluxSurface : public FESurface
{
public:
	FEHeatFluxSurface(FEMesh* pm) : FESurface(pm){}

	//! allocate storage
	void create(int n)
	{
		FESurface::create(n);
		m_FC.resize(n);
	}

	//! clone
	FEDomain* Clone()
	{
		FEHeatFluxSurface* ps = new FEHeatFluxSurface(m_pMesh);
		ps->m_FC = m_FC;
		return ps;
	}

	//! get a heat flux load BC
	FEHeatFlux& HeatFlux(int n) { return m_FC[n]; }

protected:
	vector<FEHeatFlux>	m_FC;
};
