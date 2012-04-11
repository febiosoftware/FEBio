#pragma once

#include "FECore/FEContactInterface.h"
#include "FEContactSurface.h"

//-----------------------------------------------------------------------------

class FEFacetSlidingSurface : public FEContactSurface
{
public:
	//! constructor
	FEFacetSlidingSurface(FEMesh* pm) : FEContactSurface(pm) { m_NQ.Attach(this); }

	//! initialization
	bool Init();

	//! Find element that contains the projection of x
	FEElement* FindMasterSegment(vec3d& x, vec3d& q, vec2d& r, bool& binit_nq, double tol, double srad);

	void ShallowCopy(FEFacetSlidingSurface& s);

	void Serialize(DumpFile& ar);

public:
	vector<double>				m_gap;	//!< gap function at integration points
	vector<vec3d>				m_nu;	//!< master normal at integration points
	vector<vec2d>				m_rs;	//!< natural coordinates of projection of integration point
	vector<double>				m_Lm;	//!< lagrange multipliers 
	vector<FESurfaceElement*>	m_pme;	//!< master element of projected integration point
	vector<int>					m_nei;	//!< surface element indices into arrays
	vector<double>				m_eps;	//!< penalty values for each integration point
	vector<double>				m_Ln;	//!< net contact pressure

protected:
	FENNQuery	m_NQ;	//!< used to find the nearest neighbour
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

	//! update 
	void Update(int niter);

	//! Create a shallow copy
	void ShallowCopy(FEContactInterface& ci);

	//! calculate contact forces
	void ContactForces(vector<double>& F, FENLSolver* psolver);

	//! calculate contact stiffness
	void ContactStiffness(FENLSolver* psolver);

	//! calculate contact pressures for file output
	void UpdateContactPressures();

	//! calculate Lagrangian augmentations
	bool Augment(int naug);

	//! serialize data to archive
	void Serialize(DumpFile& ar);

protected:
	//! project slave surface onto master
	void ProjectSurface(FEFacetSlidingSurface& ss, FEFacetSlidingSurface& ms);

	//! calculate auto-penalty
	void CalcAutoPenalty(FEFacetSlidingSurface& s);

public:
	double	m_epsn;			//!< normal penalty factor
	double	m_knmult;		//!< normal stiffness multiplier
	double	m_stol;			//!< search tolerance
	bool	m_btwo_pass;	//!< two-pass flag
	bool	m_bautopen;		//!< auto-penalty flag
	double	m_srad;			//!< search radius (% of model size)

	double	m_atol;			//!< aug lag tolernace
	double	m_gtol;			//!< gap tolerance
	int		m_naugmin;		//!< min nr of augmentations
	int		m_naugmax;		//!< max nr of augmentations

	double			m_mu;		//!< friction coefficient (not implemented yet)
	double			m_epsf;		//!< penalty scale factor for friction (not implementer yet)

	double	m_dxtol;		//!< penalty insertion distance

	FEFacetSlidingSurface	m_ms;	//!< master surface
	FEFacetSlidingSurface	m_ss;	//!< slave surface

	DECLARE_PARAMETER_LIST();
};
