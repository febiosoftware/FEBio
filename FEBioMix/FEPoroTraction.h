#pragma once
#include "FECore/FESurfaceLoad.h"

//-----------------------------------------------------------------------------
//! This boundary condition applies a poro-elastic normal traction on a surface
//!
class FEPoroNormalTraction : public FESurfaceLoad
{
public:
	struct LOAD
	{
		LOAD() { s[0] = s[1] = s[2] = s[3] = s[4] = s[5] = s[6] = s[7] = 1.0; }
		double	s[8];		// nodal scale factors
		int		lc;			// load curve
	};

public:
	//! constructor
	FEPoroNormalTraction(FESurface* ps, bool blinear = false, bool beffective = false) : FESurfaceLoad(ps) { m_blinear = blinear; m_beffective = beffective; }

	//! allocate storage
	void create(int n) { m_PC.resize(n); }
/*
	//! clone
	FEDomain* Clone()
	{
		FEPoroNormalTraction* ps = new FEPoroNormalTraction(m_pMesh);
		ps->m_PC = m_PC;
		return ps;
	}
*/
	//! get a pressure load BC
	LOAD& NormalTraction(int n) { return m_PC[n]; }

	//! calculate pressure stiffness
	void StiffnessMatrix(FESolver* psolver);

	//! calculate residual
	void Residual(FEGlobalVector& R);

	//! serialize data
	void Serialize(DumpFile& ar);

protected:
	//! calculate stiffness for an element
	void TractionStiffness(FESurfaceElement& el, matrix& ke, vector<double>& tn, bool effective, bool bsymm);

	//! Calculates external pressure forces
	bool TractionForce(FESurfaceElement& el, vector<double>& fe, vector<double>& tn);

	//! Calculates the linear external pressure forces (ie. non-follower forces)
	bool LinearTractionForce(FESurfaceElement& el, vector<double>& fe, vector<double>& tn);

protected:
	bool	m_blinear;		//!< linear or not (true is non-follower, false is follower)
	bool	m_beffective;	//!< effective or total normal traction

	// pressure boundary data
	vector<LOAD>	m_PC;		//!< pressure boundary cards
};
