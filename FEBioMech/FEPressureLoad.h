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
		double	s[8];		// nodal scale factors
		int		lc;			// load curve

		LOAD() { s[0] = s[1] = s[2] = s[3] = s[4] = s[5] = s[6] = s[7] = 1.0; }
	};

public:
	//! constructor
	FEPressureLoad(FEModel* pfem) : FESurfaceLoad(pfem) { m_blinear = false; }

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
	void StiffnessMatrix(FESolver* psolver);

	//! calculate residual
	void Residual(FEGlobalVector& R);

	//! serialize data
	void Serialize(DumpFile& ar);

	//! set the linear flag
	void SetLinear(bool blinear) { m_blinear = blinear; }

	//! Check if this is a linear force or not
	bool IsLinear() { return m_blinear; }

public:
	//! set an attribute of the surface load
	bool SetAttribute(const char* szatt, const char* szval);

	//! set an attribute of a surface facet
	bool SetFacetAttribute(int nface, const char* szatt, const char* szval);

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
