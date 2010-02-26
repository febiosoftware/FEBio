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
	FEPressureLoad() { s[0] = s[1] = s[2] = s[3] = 1.0; bc = 0; blinear = false; }

public:
	double	s[4];		// nodal scale factors
	int		face;		// face number
	int		lc;			// load curve
	int		bc;			// degree of freedom
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
	void PressureStiffness(FESurfaceElement& el, matrix& ke);

	//! Calculates external pressure forces
	bool PressureForce(FESurfaceElement& el, vector<double>& fe);

	//! Calculates the linear external pressure forces (ie. non-follower forces)
	bool LinearPressureForce(FESurfaceElement& el, vector<double>& fe);

protected:
	// pressure boundary data
	vector<FEPressureLoad>	m_PC;		//!< pressure boundary cards
};
