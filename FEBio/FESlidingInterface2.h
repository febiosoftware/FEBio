#pragma once
#include "FEContactInterface.h"
#include "vec2d.h"

//-----------------------------------------------------------------------------
class FEContactSurface2 : public FESurface
{
public:
	//! constructor
	FEContactSurface2(FEMesh* pm) : FESurface(pm) {}

	//! initialization
	void Init();

	//! shallow copy
	void ShallowCopy(FEContactSurface2& s);

	//! find the intersection of a ray with the surface
	FESurfaceElement* FindIntersection(vec3d r, vec3d n, double rs[2]);

public:
	vector<double>				m_gap;	//!< gap function at integration points
	vector<vec3d>				m_nu;	//!< normal at integration points
	vector<vec2d>				m_rs;	//!< natural coordinates of projection of integration point
	vector<double>				m_Lm;	//!< lagrange multipliers 
	vector<FESurfaceElement*>	m_pme;	//!< master element of projected integration point
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

	int		m_npass;	//!< nr of passes
	double	m_atol;		//!< augmentation tolerance
	double	m_eps;		//!< penalty factor
};
