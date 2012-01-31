#pragma once
#include "FECore/FESurfaceLoad.h"

//-----------------------------------------------------------------------------
//! The pressure surface is a surface domain that sustains pressure boundary
//! conditions
//!
class FEPressureLoad : public FESurfaceLoad
{
public:
	struct LOAD
	{
		double	s[4];		// nodal scale factors
		int		lc;			// load curve

		LOAD() { s[0] = s[1] = s[2] = s[3] = 1.0; }
	};

public:
	//! constructor
	FEPressureLoad(FESurface* ps, bool blinear = false) : FESurfaceLoad(ps) { m_blinear = blinear; }

	//! allocate storage
	void create(int n) { m_PC.resize(n); }

	//! clone
/*	FEDomain* Clone()
	{
		FEPressureLoad* ps = new FEPressureLoad(m_surf.GetMesh());
		ps->m_PC = m_PC;
		return ps;
	}
*/

	//! get a pressure load BC
	LOAD& PressureLoad(int n) { return m_PC[n]; }

	//! calculate pressure stiffness
	void StiffnessMatrix(FENLSolver* psolver);

	//! calculate residual
	void Residual(FENLSolver* psolver, vector<double>& R);

	//! serialize data
	void Serialize(DumpFile& ar);

	//! Check if this is a linear force or not
	bool IsLinear() { return m_blinear; }

protected:
	//! calculate stiffness for an element
	void PressureStiffness(FESurfaceElement& el, matrix& ke, vector<double>& tn);

	//! Calculates external pressure forces
	bool PressureForce(FESurfaceElement& el, vector<double>& fe, vector<double>& tn);

	//! Calculates the linear external pressure forces (ie. non-follower forces)
	bool LinearPressureForce(FESurfaceElement& el, vector<double>& fe, vector<double>& tn);

protected:
	bool			m_blinear;	//!< pressure load type (linear or nonlinear)
	vector<LOAD>	m_PC;		//!< pressure load cards
};
