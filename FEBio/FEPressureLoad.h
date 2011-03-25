#pragma once
#include "FESurface.h"
#include "FEBoundaryCondition.h"
#include "FESurfaceLoad.h"

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
		int		face;		// face number
		int		lc;			// load curve

		LOAD() { s[0] = s[1] = s[2] = s[3] = 1.0; }
	};

	// pressure load types
	enum { LINEAR, NONLINEAR };

public:
	//! constructor
	FEPressureLoad(FEMesh* pm) : m_surf(pm) { m_ntype = NONLINEAR; }

	//! allocate storage
	void create(int n)
	{
		m_surf.create(n);
		m_PC.resize(n);
	}

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
	void StiffnessMatrix(FESolidSolver* psolver);

	//! calculate residual
	void Residual(FESolidSolver* psolver, vector<double>& R);

	//! serialize data
	void Serialize(FEM& fem, DumpFile& ar);

	//! Surface data
	FESurface& Surface() { return m_surf; }

	//! Set the load type
	void SetType(int ntype) { m_ntype = ntype; }

protected:
	//! calculate stiffness for an element
	void PressureStiffness(FESurfaceElement& el, matrix& ke, vector<double>& tn);

	//! Calculates external pressure forces
	bool PressureForce(FESurfaceElement& el, vector<double>& fe, vector<double>& tn);

	//! Calculates the linear external pressure forces (ie. non-follower forces)
	bool LinearPressureForce(FESurfaceElement& el, vector<double>& fe, vector<double>& tn);

protected:
	int				m_ntype;	//!< pressure load type (linear or nonlinear)
	FESurface		m_surf;		//!< surface to which this BC is applied
	vector<LOAD>	m_PC;		//!< pressure load cards
};

//-----------------------------------------------------------------------------
//! FEConstTractionSurface is a surface that has a constant (deformation independant)
//! traction force on it.
//!

class FETractionLoad : public FEBoundaryCondition
{
public:
	FETractionLoad() { nface = -1; }

public:
	vec3d	s[4];		// nodal scale factors
	int		nface;		// face number
	int		lc;			// load curve
};

class FEConstTractionSurface : public FESurface
{
public:
	//! constructor
	FEConstTractionSurface(FEMesh* pm) : FESurface(pm) {}

	//! allocate storage
	void create(int n)
	{
		FESurface::create(n);
		m_TC.resize(n);
	}

	//! clone
	FEDomain* Clone()
	{
		FEConstTractionSurface* ps = new FEConstTractionSurface(m_pMesh);
		ps->m_TC = m_TC;
		return ps;
	}

	//! get a traction load BC
	FETractionLoad& TractionLoad(int n) { return m_TC[n]; }

	//! calculate pressure stiffness
	void StiffnessMatrix(FESolidSolver* psolver) {}

	//! calculate residual
	void Residual(FESolidSolver* psolver, vector<double>& R);

protected:
	// traction boundary data
	vector<FETractionLoad>	m_TC;	//!< traction boundary cards
};
