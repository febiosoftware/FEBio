#pragma once
#include "FECore/FEContactInterface.h"
#include "FEContactSurface.h"

//-----------------------------------------------------------------------------
class FETiedBiphasicSurface : public FEContactSurface
{
public:
	//! constructor
	FETiedBiphasicSurface(FEModel* pfem);
	
	//! initialization
	bool Init();
	
	//! shallow copy
	void ShallowCopy(FETiedBiphasicSurface& s);
	
	//! calculate the nodal normals
	void UpdateNodeNormals();
	
	void Serialize(DumpFile& ar);
	
	void SetPoroMode(bool bporo) { m_bporo = bporo; }
	
protected:
	FEModel*	m_pfem;
	
public:
	bool				m_bporo;	//!< set poro-mode
	
	vector<vec3d>				m_Gap;	//!< initial gap in reference configuration
	vector<vec3d>				m_dg;	//!< gap function at integration points
	vector<vec3d>				m_nu;	//!< normal at integration points
	vector<vec2d>				m_rs;	//!< natural coordinates of projection of integration point
	vector<vec3d>				m_Lmd;	//!< lagrange multipliers for displacements
	vector<double>				m_Lmp;	//!< lagrange multipliers for fluid pressures
	vector<FESurfaceElement*>	m_pme;	//!< master element of projected integration point
	vector<int>					m_nei;	//!< surface element indices into arrays
	vector<bool>				m_poro;	//!< surface element poro status
	
	vector<double>				m_epsn;	//!< penalty factors
	vector<double>				m_epsp;	//!< pressure penalty factors
	
	vector<vec3d>				m_nn;	//!< node normals
	
	// biphasic data
	vector<double>				m_pg;	//!< pressure "gap"
	
};

//-----------------------------------------------------------------------------
class FETiedBiphasicInterface :	public FEContactInterface
{
public:
	//! constructor
	FETiedBiphasicInterface(FEModel* pfem);
	
	//! destructor
	~FETiedBiphasicInterface();
	
	//! initialization
	bool Init();
	
	//! interface activation
	void Activate();

	//! update
	void Update(int n);
	
	//! Create a shallow copy
	void ShallowCopy(FEContactInterface& ci);
	
	//! calculate contact forces
	void ContactForces(FEGlobalVector& R);
	
	//! calculate contact stiffness
	void ContactStiffness(FENLSolver* psolver);
	
	//! calculate Lagrangian augmentations
	bool Augment(int naug);
	
	//! serialize data to archive
	void Serialize(DumpFile& ar);

	//! return the master and slave surface
	FESurface* GetMasterSurface() { return &m_ms; }
	FESurface* GetSlaveSurface () { return &m_ss; }
	
protected:
	void InitialProjection(FETiedBiphasicSurface& ss, FETiedBiphasicSurface& ms);
	void ProjectSurface(FETiedBiphasicSurface& ss, FETiedBiphasicSurface& ms);
	
	//! calculate penalty factor
	void CalcAutoPenalty(FETiedBiphasicSurface& s);
	
	void CalcAutoPressurePenalty(FETiedBiphasicSurface& s);
	double AutoPressurePenalty(FESurfaceElement& el, FETiedBiphasicSurface& s);
	
public:
	FETiedBiphasicSurface	m_ms;	//!< master surface
	FETiedBiphasicSurface	m_ss;	//!< slave surface
	
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
	
	double			m_epsn;			//!< normal penalty factor
	bool			m_bautopen;		//!< use autopenalty factor
	
	// biphasic contact parameters
	double			m_epsp;		//!< flow rate penalty
	
	DECLARE_PARAMETER_LIST();
};
