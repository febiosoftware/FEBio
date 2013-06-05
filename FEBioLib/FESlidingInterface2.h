#pragma once
#include "FECore/FEContactInterface.h"
#include "FEContactSurface.h"

//-----------------------------------------------------------------------------
class FESlidingSurface2 : public FEContactSurface
{
public:
	//! Integration point data
	class Data
	{
	public:
		Data();

	public:
		double	m_gap;	//!< gap function
		double	m_Lmd;	//!< lagrange multipliers for displacement
		double	m_Lmp;	//!< lagrange multipliers for fluid pressures
		double	m_Ln;	//!< net contact pressure
		double	m_epsn;	//!< penalty factor
		double	m_epsp;	//!< pressure penatly factor
		double	m_pg;	//!< pressure "gap" for biphasic contact
		vec3d	m_nu;	//!< normal at integration points
		vec2d	m_rs;	//!< natrual coordinates of projection
		FESurfaceElement*	m_pme;	//!< master element
	};

public:
	//! constructor
	FESlidingSurface2(FEModel* pfem);

	//! initialization
	bool Init();

	//! shallow copy
	void ShallowCopy(FESlidingSurface2& s);

	//! calculate the nodal normals
	void UpdateNodeNormals();

	void Serialize(DumpFile& ar);

	void SetPoroMode(bool bporo) { m_bporo = bporo; }

protected:
	FEModel*	m_pfem;

public:
	bool	m_bporo;	//!< set poro-mode

	vector< vector<Data> >	m_Data;	//!< integration point data
	vector<bool>		m_poro;	//!< surface element poro status
	vector<vec3d>		m_nn;	//!< node normals
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
	bool Init();

	//! interface activation
	void Activate();

	//! update
	void Update(int niter);

	//! Create a shallow copy
	void ShallowCopy(FEContactInterface& ci);

	//! calculate contact forces
	void ContactForces(FEGlobalVector& R);

	//! calculate contact stiffness
	void ContactStiffness(FENLSolver* psolver);

	//! calculate contact pressures for file output
	void UpdateContactPressures();

	//! calculate Lagrangian augmentations
	bool Augment(int naug);

	//! serialize data to archive
	void Serialize(DumpFile& ar);

	//! mark free-draining condition 
	void MarkFreeDraining();
	
	//! set free-draining condition 
	void SetFreeDraining();

	//! return the master and slave surface
	FESurface* GetMasterSurface() { return &m_ms; }
	FESurface* GetSlaveSurface () { return &m_ss; }

	
protected:
	void ProjectSurface(FESlidingSurface2& ss, FESlidingSurface2& ms, bool bupseg);

	//! calculate penalty factor
	void CalcAutoPenalty(FESlidingSurface2& s);

	void CalcAutoPressurePenalty(FESlidingSurface2& s);
	double AutoPressurePenalty(FESurfaceElement& el, FESlidingSurface2& s);

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
