#pragma once
#include "FEContactInterface.h"
#include "FEContactSurface.h"

//-----------------------------------------------------------------------------
//! Surface definition for the facet-to-facet tied interface
class FEFacetTiedSurface : public FEContactSurface
{
public:
	//! integration point data
	class Data
	{
	public:
		Data();

	public:
		vec3d	m_gap;	//!< gap function
		vec3d	m_Lm;	//!< Lagrange multiplier
		vec2d	m_rs;	//!< natural coordinates on master element
		FESurfaceElement*	m_pme;	//!< master element
	};

public:
	//! constructor
	FEFacetTiedSurface(FEModel* pfem);

	//! Initialization
	bool Init();

	//! create a shallow copy for running restarts
	void ShallowCopy(DumpStream& dmp, bool bsave);

	//! serialization for cold restarts
	void Serialize(DumpFile& ar);

public:
	vector< vector<Data> >	m_Data;	//!< integration point data
};

//-----------------------------------------------------------------------------
//! Tied contact interface with facet-to-facet integration
class FEFacet2FacetTied : public FEContactInterface
{
public:
	//! constructor
	FEFacet2FacetTied(FEModel* pfem);

	//! Initialization
	bool Init();

	//! interface activation
	void Activate();

	//! Create a shallow copy
	void ShallowCopy(DumpStream& dmp, bool bsave);

	//! serialize data to archive
	void Serialize(DumpFile& ar);

	//! return the master and slave surface
	FESurface* GetMasterSurface() { return &m_ms; }
	FESurface* GetSlaveSurface () { return &m_ss; }

	//! return integration rule class
	bool UseNodalIntegration() { return false; }

	//! build the matrix profile for use in the stiffness matrix
	void BuildMatrixProfile(FEStiffnessMatrix& K);

public:
	//! calculate contact forces
	void ContactForces(FEGlobalVector& R);

	//! calculate contact stiffness
	void ContactStiffness(FESolver* psolver);

	//! calculate Lagrangian augmentations
	bool Augment(int naug);

	//! update contact data
	void Update(int niter);

protected:

	//! projects slave nodes onto master nodes
	void ProjectSurface(FEFacetTiedSurface& ss, FEFacetTiedSurface& ms);

private:
	FEFacetTiedSurface	m_ss;	//!< slave surface
	FEFacetTiedSurface	m_ms;	//!< master surface

public:
	double		m_atol;		//!< augmentation tolerance
	double		m_eps;		//!< penalty scale factor
	double		m_stol;		//!< search tolerance
	int			m_naugmax;	//!< maximum nr of augmentations
	int			m_naugmin;	//!< minimum nr of augmentations

	DECLARE_PARAMETER_LIST();
};
