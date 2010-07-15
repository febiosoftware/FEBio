#pragma once
#include "FESurface.h"
#include "FEBoundaryCondition.h"

//-----------------------------------------------------------------------------
//! This class describes a pressure load on a surface element
//! bc = 0 for pressure loads
//! bc = 1 for heat flux

class FEPressureLoad : public FEBoundaryCondition
{
public:
	FEPressureLoad() { s[0] = s[1] = s[2] = s[3] = 1.0;
						blinear = false;}

public:
	double	s[4];		// nodal scale factors
	int		face;		// face number
	int		lc;			// load curve
	bool	blinear;	// linear or not (true is non-follower, false is follower)
};

//-----------------------------------------------------------------------------
//! The pressure surface is a surface domain that sustains pressure boundary
//! conditions
//!
class FEPressureSurface : public FESurface
{
public:
	//! constructor
	FEPressureSurface(FEMesh* pm) : FESurface(pm) {}

	//! allocate storage
	void create(int n)
	{
		FESurface::create(n);
		m_PC.resize(n);
	}

	//! clone
	FEDomain* Clone()
	{
		FEPressureSurface* ps = new FEPressureSurface(m_pMesh);
		ps->m_PC = m_PC;
		return ps;
	}

	//! get a pressure load BC
	FEPressureLoad& PressureLoad(int n) { return m_PC[n]; }

	//! calculate pressure stiffness
	void StiffnessMatrix(FESolidSolver* psolver);

	//! calculate residual
	void Residual(FESolidSolver* psolver, vector<double>& R);

	//! serialize data
	void Serialize(FEM& fem, Archive& ar);

protected:
	//! calculate stiffness for an element
	void PressureStiffness(FESurfaceElement& el, matrix& ke, double* tn);

	//! Calculates external pressure forces
	bool PressureForce(FESurfaceElement& el, vector<double>& fe, double* tn);

	//! Calculates the linear external pressure forces (ie. non-follower forces)
	bool LinearPressureForce(FESurfaceElement& el, vector<double>& fe, double* tn);

protected:
	// pressure boundary data
	vector<FEPressureLoad>	m_PC;		//!< pressure boundary cards
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
