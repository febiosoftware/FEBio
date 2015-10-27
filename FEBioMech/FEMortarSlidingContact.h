#pragma once
#include "FEContactInterface.h"
#include "FEContactSurface.h"

//-----------------------------------------------------------------------------
//! This class represents a surface used by the mortar contact interface.
class FEMortarSurface : public FEContactSurface
{
public:
	FEMortarSurface(FEMesh* pm = 0);

	//! Initializes data structures
	bool Init();

	//! update the normals
	void UpdateNormals();

public:
	vector<double>	m_p;	//!< nodal contact pressures
	vector<vec3d>	m_nu;	//!< nodal normals
	vector<vec3d>	m_gap;	//!< nodal gap function
};

//-----------------------------------------------------------------------------
class Patch
{
public:
	struct FACET
	{
		FACET(){}
		FACET(vec3d x[3]) { r[0] = x[0]; r[1] = x[1]; r[2] = x[2]; }
		vec3d	r[3];	//!< position of nodes

		// Evaluate the spatial position using iso-parametric coordinates
		vec3d Position(double rp, double sp)
		{
			return (r[0]*(1.0 - rp - sp) + r[1]*rp + r[2]*sp);
		}

		//! calculate the area of the patch
		double Area()
		{
			vec3d n = (r[1] - r[0])^(r[2] - r[0]);
			return n.norm()*0.5;
		}
	};

public:
	//! Clear the patch
	void Clear() { m_tri.clear(); }

	//! Add a facet to the patch
	void Add(vec3d x[3]) { m_tri.push_back(FACET(x)); }

	//! retrieve a facet
	FACET& Facet(int i) { return m_tri[i]; }

	//! facet count
	int Size() { return (int) m_tri.size(); }

	//! See if the facet is empty
	bool Empty() { return m_tri.empty(); }

private:
	vector<FACET>	m_tri;	//!< triangular patches
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
	//! Calculate intersection of two patches
	bool CalculateIntersection(int k, int l, Patch& patch);

	//! Update the gap values
	void UpdateNodalGaps();

	//! Update integration factors
	void UpdateMortarWeights();

private:
	bool	m_blaugon;	//!< augmented Lagrangian flag
	double	m_atol;		//!< augmented Lagrangian tolerance
	double	m_eps;		//!< penalty factor
	int		m_naugmin;	//!< minimum number of augmentations
	int		m_naugmax;	//!< maximum number of augmentations

private:
	FEMortarSurface	m_ms;	//!< mortar surface
	FEMortarSurface	m_ss;	//!< non-mortar surface

	matrix	m_n1;	//!< integration weights n1_AB
	matrix	m_n2;	//!< integration weights n2_AB

	// integration rule
	FESurfaceElementTraits*	m_pT;

	DECLARE_PARAMETER_LIST();
};
