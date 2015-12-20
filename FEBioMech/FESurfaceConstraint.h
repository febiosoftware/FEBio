#pragma once

#include "FEContactInterface.h"
#include "FEContactSurface.h"

//-----------------------------------------------------------------------------
class FESurfaceConstraintSurface : public FEContactSurface
{
public:
	//! constructor
	FESurfaceConstraintSurface(FEModel* pfem) : FEContactSurface(pfem) { m_nref = -1; }

	//! initializes data
	bool Init();

	//! shallow copy
	void ShallowCopy(DumpStream& dmp, bool bsave);

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

	//! interface activation
	void Activate();

	//! update
	void Update(int niter);

	//! shallow copy
	void ShallowCopy(DumpStream& dmp, bool bsave);

	//! calculate contact forces
	void ContactForces(FEGlobalVector& R);

	//! calculate contact stiffness
	void ContactStiffness(FESolver* psolver);

	//! calculate Lagrangian augmentations
	bool Augment(int naug);

	//! serialize data to archive
	void Serialize(DumpFile& ar);

	//! return the master and slave surface
	FESurface* GetMasterSurface() { return &m_ms; }
	FESurface* GetSlaveSurface () { return &m_ss; }

	//! return integration rule class
	bool UseNodalIntegration() { return true; }

	//! build the matrix profile for use in the stiffness matrix
	void BuildMatrixProfile(FEStiffnessMatrix& K);

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

	int	m_dofX;
	int m_dofY;
	int	m_dofZ;

	DECLARE_PARAMETER_LIST();
};
