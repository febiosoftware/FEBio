#pragma once
#include <FECore/FESurfaceLoad.h>

//-----------------------------------------------------------------------------
//! Surface that sustains a convective heat flux boundary condition
//!
class FEConvectiveHeatFlux : public FESurfaceLoad
{
public:
	struct LOAD
	{
		double	hc;			// heat transfer coefficient
		double	s[8];		// nodal scale factors
		int		lc;			// load curve describing ambient temperature
		LOAD()  { s[0] = s[1] = s[2] = s[3] = s[4] = s[5] = s[6] = s[7] = 1.0; lc = -1; }
	};

public:
	FEConvectiveHeatFlux(FESurface* ps) : FESurfaceLoad(ps){}

	//! allocate storage
	void create(int n) { m_FC.resize(n); }

	//! clone
/*	FEDomain* Clone()
	{
		FEConvectiveHeatFlux* ps = new FEConvectiveHeatFlux(m_pMesh);
		ps->m_FC = m_FC;
		return ps;
	}
*/

	//! get a heat flux load BC
	LOAD& HeatFlux(int n) { return m_FC[n]; }

	//! stiffness matrix
	void StiffnessMatrix(FENLSolver* psolver);
	
	//! residual
	void Residual(FEGlobalVector& R);

	//! serialization
	void Serialize(DumpFile& ar) {}

protected:
	void ElementStiffness(FESurfaceElement& el, matrix& ke, double hc);

protected:
	vector<LOAD>	m_FC;
};
