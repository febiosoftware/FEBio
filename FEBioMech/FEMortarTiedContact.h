#pragma once
#include "FEMortarInterface.h"
#include "FEMortarContactSurface.h"

//-----------------------------------------------------------------------------
//! This class represents a surface used by the mortar contact interface.
class FEMortarTiedSurface : public FEMortarContactSurface
{
public:
	FEMortarTiedSurface(FEMesh* pm = 0);

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
	FESurface* GetMasterSurface() { return &m_ms; }
	FESurface* GetSlaveSurface () { return &m_ss; }

public:
	//! temporary construct to determine if contact interface uses nodal integration rule (or facet)
	bool UseNodalIntegration() { return false; }

	//! interface activation
	void Activate();

	//! one-time initialization
	bool Init();

	//! calculate contact forces
	void ContactForces(FEGlobalVector& R);

	//! calculate contact stiffness
	void ContactStiffness(FESolver* psolver);

	//! calculate Lagrangian augmentations
	bool Augment(int naug);

	//! serialize data to archive
	void Serialize(DumpFile& ar);

	//! build the matrix profile for use in the stiffness matrix
	void BuildMatrixProfile(FEStiffnessMatrix& K);

	//! update interface data
	void Update(int niter);

	//! shallow copy
	void ShallowCopy(DumpStream& dmp, bool bsave);

private:
	double	m_atol;		//!< augmented Lagrangian tolerance
	double	m_eps;		//!< penalty factor
	int		m_naugmin;	//!< minimum number of augmentations
	int		m_naugmax;	//!< maximum number of augmentations

private:
	FEMortarTiedSurface	m_ms;	//!< mortar surface
	FEMortarTiedSurface	m_ss;	//!< non-mortar surface

	DECLARE_PARAMETER_LIST();
};
