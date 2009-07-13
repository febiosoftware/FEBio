#pragma once

#include "FEContactInterface.h"
#include "FESurface.h"
#include "vec2d.h"

//-----------------------------------------------------------------------------

class FEFacetSlidingSurface : public FESurface
{
public:
	//! constructor
	FEFacetSlidingSurface(FEMesh* pm) : FESurface(pm) { m_NQ.Attach(this); }

	//! initialization
	void Init();

	//! Find element that contains the projection of x
	FEElement* FindMasterSegment(vec3d& x, vec3d& q, vec2d& r, bool& binit_nq, double tol);

	void ShallowCopy(FEFacetSlidingSurface& s);

public:
	vector<double>				m_gap;	//!< gap function at integration points
	vector<vec3d>				m_nu;	//!< master normal at integration points
	vector<vec2d>				m_rs;	//!< natural coordinates of projection of integration point
	vector<double>				m_Lm;	//!< lagrange multipliers 
	vector<FESurfaceElement*>	m_pme;	//!< master element of projected integration point

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
	FEFacet2FacetSliding(FEM* pfem);

	//! initialization routine
	void Init();

	//! update 
	void Update();

	//! Create a shallow copy
	void ShallowCopy(FEContactInterface& ci);

	//! calculate contact forces
	void ContactForces(vector<double>& F);

	//! calculate contact stiffness
	void ContactStiffness();

	//! calculate Lagrangian augmentations
	bool Augment(int naug);

	//! serialize data to archive
	void Serialize(Archive& ar);

protected:
	//! project slave surface onto master
	void ProjectSurface(FEFacetSlidingSurface& ss, FEFacetSlidingSurface& ms);

public:
	double	m_epsn;		//!< normal penalty factor
	double	m_knmult;	//!< normal stiffness multiplier
	double	m_stol;		//!< search tolerance
	int		m_npass;	//!< nr of passes

	double	m_atol;		//!< aug lag tolernace
	double	m_gtol;		//!< gap tolerance
	int		m_naugmin;	//!< min nr of augmentations
	int		m_naugmax;	//!< max nr of augmentations

	FEFacetSlidingSurface	m_ms;	//!< master surface
	FEFacetSlidingSurface	m_ss;	//!< slave surface
};
