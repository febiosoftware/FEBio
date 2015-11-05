#pragma once
#include "FEContactInterface.h"
#include "FEContactSurface.h"
#include "FECore/mortar.h"

//-----------------------------------------------------------------------------
//! This class represents a surface used by the mortar contact interface.
class FEMortarSlidingSurface : public FEContactSurface
{
public:
	FEMortarSlidingSurface(FEMesh* pm = 0);

	//! Initializes data structures
	bool Init();

	//! update the normals
	void UpdateNormals(bool binit);

	//! update nodal areas
	void UpdateNodalAreas();

public:
	vector<double>	m_p;		//!< nodal contact pressures
	vector<double>	m_L;		//!< Lagrange multipliers
	vector<vec3d>	m_nu;		//!< nodal normals
	vector<double>	m_norm0;	//!< initial (inverse) normal lenghts
	vector<vec3d>	m_gap;		//!< nodal gap function
	vector<double>	m_A;		//!< nodal areas
};

//-----------------------------------------------------------------------------
//! This class implements a mortar contact formulation for frictionless, sliding contact
class FEMortarSlidingContact : public FEContactInterface
{
public:
	//! constructor
	FEMortarSlidingContact(FEModel* pfem);

	//! destructor
	~FEMortarSlidingContact();

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

protected:
	//! Update the gap values
	void UpdateNodalGaps();

	//! Update integration factors
	void UpdateMortarWeights();

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

	matrix	m_n1;	//!< integration weights n1_AB
	matrix	m_n2;	//!< integration weights n2_AB

	// integration rule
	FESurfaceElementTraits*	m_pT;

	DECLARE_PARAMETER_LIST();
};
