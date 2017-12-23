#pragma once
#include "FEMortarInterface.h"
#include "FEMortarContactSurface.h"

//-----------------------------------------------------------------------------
//! This class represents a surface used by the mortar contact interface.
class FEMortarTiedSurface : public FEMortarContactSurface
{
public:
	FEMortarTiedSurface(FEModel* pfem);

	//! Initializes data structures
	bool Init();

public:
	vector<vec3d>	m_L;		//!< Lagrange multipliers
};

//-----------------------------------------------------------------------------
//! This class implements a mortar based tied interface.
class FEMortarTiedContact : public FEMortarInterface
{
public:
	//! constructor
	FEMortarTiedContact(FEModel* pfem);

	//! return the master and slave surface
	FESurface* GetMasterSurface() override { return &m_ms; }
	FESurface* GetSlaveSurface () override { return &m_ss; }

public:
	//! temporary construct to determine if contact interface uses nodal integration rule (or facet)
	bool UseNodalIntegration() override { return false; }

	//! interface activation
	void Activate() override;

	//! one-time initialization
	bool Init() override;

	//! calculate contact forces
	void Residual(FEGlobalVector& R, const FETimeInfo& tp) override;

	//! calculate contact stiffness
	void StiffnessMatrix(FESolver* psolver, const FETimeInfo& tp) override;

	//! calculate Lagrangian augmentations
	bool Augment(int naug, const FETimeInfo& tp) override;

	//! serialize data to archive
	void Serialize(DumpStream& ar) override;

	//! build the matrix profile for use in the stiffness matrix
	void BuildMatrixProfile(FEGlobalMatrix& K) override;

	//! update interface data
	void Update(int niter, const FETimeInfo& tp) override;

private:
	double	m_atol;		//!< augmented Lagrangian tolerance
	double	m_eps;		//!< penalty factor
	int		m_naugmin;	//!< minimum number of augmentations
	int		m_naugmax;	//!< maximum number of augmentations

private:
	FEMortarTiedSurface	m_ms;	//!< mortar surface
	FEMortarTiedSurface	m_ss;	//!< non-mortar surface

	int		m_dofX;
	int		m_dofY;
	int		m_dofZ;

	DECLARE_PARAMETER_LIST();
};
