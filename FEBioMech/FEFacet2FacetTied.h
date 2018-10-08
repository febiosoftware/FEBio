#pragma once
#include "FEContactInterface.h"
#include "FEContactSurface.h"

//-----------------------------------------------------------------------------
//! Surface definition for the facet-to-facet tied interface
class FEFacetTiedSurface : public FEContactSurface
{
public:
	//! integration point data
	class Data : public FEContactMaterialPoint
	{
	public:
		Data();

	public:
		vec3d	m_vgap;	//!< gap function
		vec3d	m_Lm;	//!< Lagrange multiplier
		vec2d	m_rs;	//!< natural coordinates on master element
		FESurfaceElement*	m_pme;	//!< master element
	};

public:
	//! constructor
	FEFacetTiedSurface(FEModel* pfem);

	//! Initialization
	bool Init();

	//! serialization for cold restarts
	void Serialize(DumpStream& ar);

	//! create material point data
	FEMaterialPoint* CreateMaterialPoint() override;
};

//-----------------------------------------------------------------------------
//! Tied contact interface with facet-to-facet integration
class FEFacet2FacetTied : public FEContactInterface
{
public:
	//! constructor
	FEFacet2FacetTied(FEModel* pfem);

	//! Initialization
	bool Init() override;

	//! interface activation
	void Activate() override;

	//! serialize data to archive
	void Serialize(DumpStream& ar) override;

	//! return the master and slave surface
	FESurface* GetMasterSurface() override { return &m_ms; }
	FESurface* GetSlaveSurface () override { return &m_ss; }

	//! return integration rule class
	bool UseNodalIntegration() override { return false; }

	//! build the matrix profile for use in the stiffness matrix
	void BuildMatrixProfile(FEGlobalMatrix& K) override;

public:
	//! calculate contact forces
	void Residual(FEGlobalVector& R, const FETimeInfo& tp) override;

	//! calculate contact stiffness
	void StiffnessMatrix(FESolver* psolver, const FETimeInfo& tp) override;

	//! calculate Lagrangian augmentations
	bool Augment(int naug, const FETimeInfo& tp) override;

	//! update contact data
	void Update(int niter, const FETimeInfo& tp) override;

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

	DECLARE_FECORE_CLASS();
};
