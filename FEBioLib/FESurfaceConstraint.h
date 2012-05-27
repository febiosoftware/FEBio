#pragma once

#include "FECore/FEContactInterface.h"
#include "FEContactSurface.h"

//-----------------------------------------------------------------------------
class FESurfaceConstraintSurface : public FEContactSurface
{
public:
	//! constructor
	FESurfaceConstraintSurface(FEMesh* pm = 0) : FEContactSurface(pm) { m_nref = -1; }

	//! initializes data
	bool Init();

	//! shallow copy
	void ShallowCopy(FESurfaceConstraintSurface& s)
	{
		m_Lm = s.m_Lm;
		m_gap = s.m_gap;
	}

	//! update surface data
	void Update();

	//! calculates the center of mass of the surface
	vec3d CenterOfMass();

	void Serialize(DumpFile& ar);

public:
	vector<vec3d>				m_gap;	//!< gap function at nodes
	vector<FESurfaceElement*>	m_pme;	//!< master element a slave node penetrates
	vector<vec2d>				m_rs;	//!< natural coordinates of slave projection on master element
	vector<vec3d>				m_Lm;	//!< Lagrange multipliers

	int		m_nref;	//!< reference node
};

//-----------------------------------------------------------------------------

class FESurfaceConstraint : public FEContactInterface
{
public:
	//! constructor
	FESurfaceConstraint(FEModel* pfem);

	//! destructor
	virtual ~FESurfaceConstraint(void) {}

	//! initialization
	bool Init();

	//! update
	void Update(int niter);

	//! shallow copy
	void ShallowCopy(FEContactInterface& ci);

	//! calculate contact forces
	void ContactForces(vector<double>& F, FENLSolver* psolver);

	//! calculate contact stiffness
	void ContactStiffness(FENLSolver* psolver);

	//! calculate Lagrangian augmentations
	bool Augment(int naug);

	//! serialize data to archive
	void Serialize(DumpFile& ar);

protected:
	void ProjectSurface(FESurfaceConstraintSurface& ss, FESurfaceConstraintSurface& ms, bool bmove);

public:
	FESurfaceConstraintSurface	m_ss;	//!< slave surface
	FESurfaceConstraintSurface	m_ms;	//!< master surface

	double	m_atol;			//!< augmentation tolerance
	double	m_eps;			//!< penalty scale factor
	double	m_stol;			//!< search tolerance
	double	m_srad;			//!< search radius (%)
	bool	m_btwo_pass;	//!< nr of passes

	DECLARE_PARAMETER_LIST();
};
