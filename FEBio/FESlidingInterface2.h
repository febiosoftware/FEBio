#pragma once
#include "FEContactInterface.h"
#include "FEContactSurface.h"
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

	//! calculate the nodal normals
	void UpdateNodeNormals();

protected:
	FEM*	m_pfem;

public:
	vector<double>				m_gap;	//!< gap function at integration points
	vector<vec3d>				m_nu;	//!< normal at integration points
	vector<vec2d>				m_rs;	//!< natural coordinates of projection of integration point
	vector<double>				m_Lmd;	//!< lagrange multipliers for displacements
	vector<double>				m_Lmp;	//!< lagrange multipliers for fluid pressures
	vector<FESurfaceElement*>	m_pme;	//!< master element of projected integration point
	vector<int>					m_nei;	//!< surface element indices into arrays

	vector<vec3d>		m_nn;	//!< node normals

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
	~FESlidingInterface2();

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
	double AutoPressurePenalty(FESurface& ss, FESurface& ms);

public:
	FEContactSurface2	m_ms;	//!< master surface
	FEContactSurface2	m_ss;	//!< slave surface

	int				m_knmult;	//!< higher order stiffness multiplier
	int				m_npass;	//!< nr of passes
	double			m_atol;		//!< augmentation tolerance
	double			m_gtol;		//!< gap tolerance
	double			m_ptol;		//!< pressure gap tolerance
	double			m_stol;		//!< search tolerance
	bool			m_bsymm;	//!< use symmetric stiffness components only
	double			m_srad;		//!< contact search radius
	int				m_naugmax;	//!< maximum nr of augmentations
	int				m_naugmin;	//!< minimum nr of augmentations

	double			m_eps;		//!< penalty factor
	bool			m_bautopen;	//!< use autopenalty factor

	bool	m_bdebug;		// debug flag
	char	m_szdebug[256];	// debug file name
	FILE*	m_fp;			// debug file
	
	// bihpasic contact parameters
	double	m_epsp;		//!< flow rate penalty
};
