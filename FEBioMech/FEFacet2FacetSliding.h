#pragma once

#include "FEContactInterface.h"
#include "FEContactSurface.h"

//-----------------------------------------------------------------------------
//! Contact surface for facet-to-facet sliding interfaces
class FEFacetSlidingSurface : public FEContactSurface
{
public:
	//! Integration point data
	class Data
	{
	public:
		Data();

	public:
		double	m_gap;	//!< gap function at integration points
		double	m_Lm;	//!< Lagrange multipliers
		double	m_eps;	//!< penalty value at integration point
		double	m_Ln;	//!< net contact pressure
		vec3d	m_nu;	//!< master normal at integration points
		vec2d	m_rs;	//!< natural coordinates of projection of integration point
		FESurfaceElement*	m_pme;	//!< master element of projection
	};

public:
	//! constructor
	FEFacetSlidingSurface(FEMesh* pm) : FEContactSurface(pm) {}

	//! initialization
	bool Init();

	//! create a shallow copy for running restarts
	void ShallowCopy(DumpStream& dmp, bool bsave);

	//! evaluate net contact force
	vec3d GetContactForce();

	//! serialize data for (cold) restart
	void Serialize(DumpFile& ar);

public:
	void GetNodalContactGap     (int nface, double* gn);
	void GetNodalContactPressure(int nface, double* pn);
	void GetNodalContactTraction(int nface, vec3d* tn);

public:
	vector< vector<Data> >		m_Data;	//!< integration point data
};

//-----------------------------------------------------------------------------
//! Sliding interface with facet to facet integration

//! This class is similar to the sliding interface except that it uses
//! a Gaussian quadrature rule in stead of a nodal integration rule
//! as its base class does.
//
class FEFacet2FacetSliding : public FEContactInterface
{
public:
	//! constructor
	FEFacet2FacetSliding(FEModel* pfem);

	//! initialization routine
	bool Init();

	//! interface activation
	void Activate();

	//! update 
	void Update(int niter);

	//! Create a shallow copy
	void ShallowCopy(DumpStream& dmp, bool bsave);

	//! calculate contact forces
	void ContactForces(FEGlobalVector& R);

	//! calculate contact stiffness
	void ContactStiffness(FESolver* psolver);

	//! calculate contact pressures for file output
	void UpdateContactPressures();

	//! calculate Lagrangian augmentations
	bool Augment(int naug);

	//! serialize data to archive
	void Serialize(DumpFile& ar);

	//! return the master and slave surface
	FESurface* GetMasterSurface() { return &m_ms; }
	FESurface* GetSlaveSurface () { return &m_ss; }

	//! return integration rule class
	bool UseNodalIntegration() { return false; }

	//! build the matrix profile for use in the stiffness matrix
	void BuildMatrixProfile(FEStiffnessMatrix& K);

protected:
	//! project slave surface onto master
	void ProjectSurface(FEFacetSlidingSurface& ss, FEFacetSlidingSurface& ms, bool bsegup);

	//! calculate auto-penalty
	void CalcAutoPenalty(FEFacetSlidingSurface& s);

public:
	double	m_epsn;			//!< normal penalty factor
	double	m_knmult;		//!< normal stiffness multiplier
	double	m_stol;			//!< search tolerance
	bool	m_btwo_pass;	//!< two-pass flag
	bool	m_bautopen;		//!< auto-penalty flag
	double	m_srad;			//!< search radius (% of model size)
	int		m_nsegup;		//!< segment update parameter

	double	m_atol;			//!< aug lag tolernace
	double	m_gtol;			//!< gap tolerance
	int		m_naugmin;		//!< min nr of augmentations
	int		m_naugmax;		//!< max nr of augmentations

	double			m_mu;		//!< friction coefficient (not implemented yet)
	double			m_epsf;		//!< penalty scale factor for friction (not implementer yet)

	double	m_dxtol;		//!< penalty insertion distance

	FEFacetSlidingSurface	m_ms;	//!< master surface
	FEFacetSlidingSurface	m_ss;	//!< slave surface

private:
	bool	m_bfirst;
	double	m_normg0;

public:
	DECLARE_PARAMETER_LIST();
};
