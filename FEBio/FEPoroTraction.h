#pragma once
#include "FESurface.h"
#include "FEBoundaryCondition.h"

//-----------------------------------------------------------------------------
//! This class describes a normal traction on a porous surface element

class FEPoroNormalTraction : public FEBoundaryCondition
{
public:
	FEPoroNormalTraction() { s[0] = s[1] = s[2] = s[3] = 1.0; 
						blinear = false; effective = false;}

public:
	double	s[4];		// nodal scale factors
	int		face;		// face number
	int		lc;			// load curve
	bool	blinear;	// linear or not (true is non-follower, false is follower)
	bool	effective;	// effective or total normal traction
};

//-----------------------------------------------------------------------------
//! The pressure surface is a surface domain that sustains pressure boundary
//! conditions
//!
class FEPoroTractionSurface : public FESurface
{
public:
	//! constructor
	FEPoroTractionSurface(FEMesh* pm) : FESurface(pm) {}

	//! allocate storage
	void create(int n)
	{
		FESurface::create(n);
		m_PC.resize(n);
	}

	//! clone
	FEDomain* Clone()
	{
		FEPoroTractionSurface* ps = new FEPoroTractionSurface(m_pMesh);
		ps->m_PC = m_PC;
		return ps;
	}

	//! get a pressure load BC
	FEPoroNormalTraction& NormalTraction(int n) { return m_PC[n]; }

	//! calculate pressure stiffness
	void StiffnessMatrix(FESolidSolver* psolver);

	//! calculate residual
	void Residual(FESolidSolver* psolver, vector<double>& R);

	//! serialize data
	void Serialize(FEM& fem, DumpFile& ar);

protected:
	//! calculate stiffness for an element
	void TractionStiffness(FESurfaceElement& el, matrix& ke, vector<double>& tn, bool effective);

	//! Calculates external pressure forces
	bool TractionForce(FESurfaceElement& el, vector<double>& fe, vector<double>& tn);

	//! Calculates the linear external pressure forces (ie. non-follower forces)
	bool LinearTractionForce(FESurfaceElement& el, vector<double>& fe, vector<double>& tn);

protected:
	// pressure boundary data
	vector<FEPoroNormalTraction>	m_PC;		//!< pressure boundary cards
};

