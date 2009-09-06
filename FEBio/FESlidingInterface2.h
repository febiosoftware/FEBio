#pragma once
#include "FEContactInterface.h"
#include "vec2d.h"

//-----------------------------------------------------------------------------
class FEContactSurface2 : public FESurface
{
public:
	//! constructor
	FEContactSurface2(FEM* pfem);

	//! initialization
	void Init();

	//! shallow copy
	void ShallowCopy(FEContactSurface2& s);

	//! find the intersection of a ray with the surface
	FESurfaceElement* FindIntersection(vec3d r, vec3d n, double rs[2], double eps);

public:
	bool Intersect(FESurfaceElement& el, vec3d r, vec3d n, double rs[2], double& g, double eps);
	bool IntersectTri(vec3d* y, vec3d r, vec3d n, double rs[2], double& g, double eps);
	bool IntersectQuad(vec3d* y, vec3d r, vec3d n, double rs[2], double& g, double eps);

protected:
	FEM*	m_pfem;

public:
	vector<double>				m_gap;	//!< gap function at integration points
	vector<vec3d>				m_nu;	//!< normal at integration points
	vector<vec2d>				m_rs;	//!< natural coordinates of projection of integration point
	vector<double>				m_Lm;	//!< lagrange multipliers 
	vector<FESurfaceElement*>	m_pme;	//!< master element of projected integration point

	// biphasic data
	vector<double>				m_pg;	//!< pressure "gap"
};

//-----------------------------------------------------------------------------
class FESlidingInterface2 :	public FEContactInterface
{
public:
	//! constructor
	FESlidingInterface2(FEM* pfem);

	//! destructor
	~FESlidingInterface2(void) {}

	//! initialization
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
	void ProjectSurface(FEContactSurface2& ss, FEContactSurface2& ms);

public:
	FEContactSurface2	m_ms;	//!< master surface
	FEContactSurface2	m_ss;	//!< slave surface

	int		m_knmult;	//!< higher order stiffness multiplier
	int		m_npass;	//!< nr of passes
	double	m_atol;		//!< augmentation tolerance
	double	m_eps;		//!< penalty factor
	double	m_stol;		//!< search tolerance
	bool	m_bsymm;	//!< use symmetric stiffness components only

	// bihpasic contact parameters
	double	m_epsp;		//!< flow rate penalty
};
