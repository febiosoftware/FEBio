#pragma once
#include "FEMortarInterface.h"
#include "FEMortarContactSurface.h"

//-----------------------------------------------------------------------------
//! This class represents a surface used by the mortar contact interface.
class FEMortarSlidingSurface : public FEMortarContactSurface
{
public:
	FEMortarSlidingSurface(FEModel* pfem);

	//! Initializes data structures
	bool Init();

	//! update the normals
	void UpdateNormals(bool binit);

public:
	vector<double>	m_p;		//!< nodal contact pressures
	vector<double>	m_L;		//!< Lagrange multipliers
	vector<vec3d>	m_nu;		//!< nodal normals
	vector<double>	m_norm0;	//!< initial (inverse) normal lenghts
};

//-----------------------------------------------------------------------------
//! This class implements a mortar contact formulation for frictionless, sliding contact
class FEMortarSlidingContact : public FEMortarInterface
{
public:
	//! constructor
	FEMortarSlidingContact(FEModel* pfem);

	//! destructor
	~FEMortarSlidingContact();

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

protected:
	// contact stiffness contributions
	void ContactGapStiffness(FESolver* psolver);
	void ContactNormalStiffness(FESolver* psolver);

private:
	double	m_atol;		//!< augmented Lagrangian tolerance
	double	m_eps;		//!< penalty factor
	int		m_naugmin;	//!< minimum number of augmentations
	int		m_naugmax;	//!< maximum number of augmentations

private:
	FEMortarSlidingSurface	m_ms;	//!< mortar surface
	FEMortarSlidingSurface	m_ss;	//!< non-mortar surface

	int		m_dofX;
	int		m_dofY;
	int		m_dofZ;

	DECLARE_FECORE_CLASS();
};
