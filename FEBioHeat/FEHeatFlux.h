#pragma once
#include <FECore/FESurfaceLoad.h>

//-----------------------------------------------------------------------------
//! Surface that sustains a heat flux boundary condition
//!
class FEHeatFlux : public FESurfaceLoad
{
public:
	struct LOAD
	{
		double	s[8];		// nodal scale factors
		int		lc;			// load curve
		LOAD()  { s[0] = s[1] = s[2] = s[3] = s[4] = s[5] = s[6] = s[7] = 1.0; lc = -1; }
	};

public:
	FEHeatFlux(FESurface* ps) : FESurfaceLoad(ps){}

	//! allocate storage
	void create(int n) { m_FC.resize(n); }

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
	void StiffnessMatrix(FENLSolver* psolver) {}
	
	//! residual
	void Residual(FEGlobalVector& R);

	//! serialization
	void Serialize(DumpFile& ar);

protected:
	vector<LOAD>	m_FC;
};
