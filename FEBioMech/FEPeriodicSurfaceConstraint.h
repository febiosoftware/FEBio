#pragma once

#include "FEContactInterface.h"
#include "FEContactSurface.h"

//-----------------------------------------------------------------------------
class FEBIOMECH_API FEPeriodicSurfaceConstraintSurface : public FEContactSurface
{
public:
	//! constructor
	FEPeriodicSurfaceConstraintSurface(FEModel* pfem) : FEContactSurface(pfem) { m_nref = -1; }

	//! initializes data
	bool Init();

	//! calculates the center of mass of the surface
	vec3d CenterOfMass();

	void Serialize(DumpStream& ar);

public:
	vector<vec3d>				m_gap;	//!< gap function at nodes
	vector<FESurfaceElement*>	m_pme;	//!< master element a slave node penetrates
	vector<vec2d>				m_rs;	//!< natural coordinates of slave projection on master element
	vector<vec3d>				m_Lm;	//!< Lagrange multipliers

	int		m_nref;	//!< reference node
};

//-----------------------------------------------------------------------------

class FEBIOMECH_API FEPeriodicSurfaceConstraint : public FEContactInterface
{
public:
	//! constructor
	FEPeriodicSurfaceConstraint(FEModel* pfem);

	//! destructor
	virtual ~FEPeriodicSurfaceConstraint(void) {}

	//! initialization
	bool Init() override;

	//! interface activation
	void Activate() override;

	//! serialize data to archive
	void Serialize(DumpStream& ar) override;

	//! return the master and slave surface
	FESurface* GetMasterSurface() override { return &m_ms; }
	FESurface* GetSlaveSurface() override { return &m_ss; }

	//! return integration rule class
	bool UseNodalIntegration() override { return true; }

	//! build the matrix profile for use in the stiffness matrix
	void BuildMatrixProfile(FEGlobalMatrix& K) override;

public:
	//! calculate contact forces
	void Residual(FEGlobalVector& R, const FETimeInfo& tp) override;

	//! calculate contact stiffness
	void StiffnessMatrix(FESolver* psolver, const FETimeInfo& tp) override;

	//! calculate Lagrangian augmentations
	bool Augment(int naug, const FETimeInfo& tp) override;

	//! update
	void Update() override;

protected:
	void ProjectSurface(FEPeriodicSurfaceConstraintSurface& ss, FEPeriodicSurfaceConstraintSurface& ms, bool bmove);

public:
	FEPeriodicSurfaceConstraintSurface	m_ss;	//!< slave surface
	FEPeriodicSurfaceConstraintSurface	m_ms;	//!< master surface

	double	m_atol;			//!< augmentation tolerance
	double	m_eps;			//!< penalty scale factor
	double	m_stol;			//!< search tolerance
	double	m_srad;			//!< search radius (%)
	bool	m_btwo_pass;	//!< nr of passes

	int	m_dofX;
	int m_dofY;
	int	m_dofZ;

	DECLARE_FECORE_CLASS();
};
