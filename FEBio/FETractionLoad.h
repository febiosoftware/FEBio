#pragma once
#include "FESurfaceLoad.h"

//-----------------------------------------------------------------------------
//! FETractionLoad is a surface that has a constant (deformation independant)
//! traction force on it.
//!
class FETractionLoad : public FESurfaceLoad
{
public:
	struct LOAD
	{
		vec3d	s[4];		// nodal scale factors
		int		lc;			// load curve
	};

public:
	//! constructor
	FETractionLoad(FEMesh* pm) : FESurfaceLoad(pm) {}

	//! allocate storage
	void create(int n)
	{
		m_surf.create(n);
		m_TC.resize(n);
	}

	//! clone
/*	FEDomain* Clone()
	{
		FETractionLoad* ps = new FETractionLoad(m_surf.GetMesh());
		ps->m_TC = m_TC;
		return ps;
	}
*/
	//! get a traction load BC
	LOAD& TractionLoad(int n) { return m_TC[n]; }

	//! calculate pressure stiffness
	void StiffnessMatrix(FESolver* psolver) {}

	//! calculate residual
	void Residual(FESolver* psolver, vector<double>& R);

protected:
	vector<LOAD>	m_TC;		//!< traction boundary cards
};
