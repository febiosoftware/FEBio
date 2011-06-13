#pragma once
#include "FECore/FEContactInterface.h"
#include "FEBioLib/FEContactSurface.h"

//-----------------------------------------------------------------------------
class FESlidingSurface2 : public FEContactSurface
{
public:
	//! constructor
	FESlidingSurface2(FEModel* pfem);

	//! initialization
	void Init();

	//! shallow copy
	void ShallowCopy(FESlidingSurface2& s);

	//! calculate the nodal normals
	void UpdateNodeNormals();

	void Serialize(DumpFile& ar);

	void SetPoroMode(bool bporo) { m_bporo = bporo; }

protected:
	FEModel*	m_pfem;

public:
	bool				m_bporo;	//!< set poro-mode

	vector<double>				m_gap;	//!< gap function at integration points
	vector<vec3d>				m_nu;	//!< normal at integration points
	vector<vec2d>				m_rs;	//!< natural coordinates of projection of integration point
	vector<double>				m_Lmd;	//!< lagrange multipliers for displacements
	vector<double>				m_Lmp;	//!< lagrange multipliers for fluid pressures
	vector<FESurfaceElement*>	m_pme;	//!< master element of projected integration point
	vector<int>					m_nei;	//!< surface element indices into arrays

	vector<double>	m_epsn;	//!< penalty factors
	vector<double>	m_epsp;	//!< pressure penalty factors

	vector<vec3d>		m_nn;	//!< node normals

	// biphasic data
	vector<double>				m_pg;	//!< pressure "gap"
};

//-----------------------------------------------------------------------------
class FESlidingInterface2 :	public FEContactInterface
{
public:
	//! constructor
	FESlidingInterface2(FEModel* pfem);

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
	void Serialize(DumpFile& ar);

	//! mark free-draining condition 
	void MarkFreeDraining();
	
	//! set free-draining condition 
	void SetFreeDraining();
	
protected:
	void ProjectSurface(FESlidingSurface2& ss, FESlidingSurface2& ms, bool bupseg);

	//! calculate penalty factor
	void CalcAutoPenalty(FESlidingSurface2& s);

	void CalcAutoPressurePenalty(FESlidingSurface2& s);
	double AutoPressurePenalty(FESurfaceElement& el, FESlidingSurface2& s);

	bool PoroStatus(FEMesh& m, FESurfaceElement& el);

public:
	FESlidingSurface2	m_ms;	//!< master surface
	FESlidingSurface2	m_ss;	//!< slave surface

	int				m_knmult;		//!< higher order stiffness multiplier
	bool			m_btwo_pass;	//!< two-pass flag
	double			m_atol;			//!< augmentation tolerance
	double			m_gtol;			//!< gap tolerance
	double			m_ptol;			//!< pressure gap tolerance
	double			m_stol;			//!< search tolerance
	bool			m_bsymm;		//!< use symmetric stiffness components only
	double			m_srad;			//!< contact search radius
	int				m_naugmax;		//!< maximum nr of augmentations
	int				m_naugmin;		//!< minimum nr of augmentations
	int				m_nsegup;		//!< segment update parameter

	double			m_epsn;		//!< normal penalty factor
	bool			m_bautopen;	//!< use autopenalty factor

	// biphasic contact parameters
	double	m_epsp;		//!< flow rate penalty

	DECLARE_PARAMETER_LIST();
};
