#pragma once
#include "FECore/FEContactInterface.h"
#include "FEContactSurface.h"

// Elastic sliding contact, reducing the algorithm of biphasic sliding contact
// (FESlidingInterface2) to elastic case.  The algorithm derives from Bonet
// & Wood's treatment of surface pressures

//-----------------------------------------------------------------------------
class FESlidingSurfaceBW : public FEContactSurface
{
public:
	//! constructor
	FESlidingSurfaceBW(FEModel* pfem);
	
	//! initialization
	bool Init();
	
	//! shallow copy
	void ShallowCopy(FESlidingSurfaceBW& s);
	
	//! calculate the nodal normals
	void UpdateNodeNormals();
	
	void Serialize(DumpFile& ar);

	//! evaluate net contact force
	vec3d NetContactForce();
	
protected:
	FEModel*	m_pfem;
	
public:
	vector<double>				m_gap;	//!< gap function at integration points
	vector<vec3d>				m_nu;	//!< normal at integration points
	vector<vec2d>				m_rs;	//!< natural coordinates of projection of integration point
	vector<double>				m_Lmd;	//!< lagrange multipliers for displacements
	vector<FESurfaceElement*>	m_pme;	//!< master element of projected integration point
	vector<int>					m_nei;	//!< surface element indices into arrays
	vector<double>				m_Ln;	//!< net contact pressure
	
	vector<double>	m_epsn;	//!< penalty factors
	
	vector<vec3d>		m_nn;	//!< node normals
	
};

//-----------------------------------------------------------------------------
class FESlidingInterfaceBW : public FEContactInterface
{
public:
	//! constructor
	FESlidingInterfaceBW(FEModel* pfem);
	
	//! destructor
	~FESlidingInterfaceBW();
	
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

	//! return the master and slave surface
	FESurface* GetMasterSurface() { return &m_ms; }
	FESurface* GetSlaveSurface () { return &m_ss; }

	
protected:
	void ProjectSurface(FESlidingSurfaceBW& ss, FESlidingSurfaceBW& ms, bool bupseg);
	
	//! calculate penalty factor
	void CalcAutoPenalty(FESlidingSurfaceBW& s);
	
public:
	FESlidingSurfaceBW	m_ms;	//!< master surface
	FESlidingSurfaceBW	m_ss;	//!< slave surface
	
	int				m_knmult;		//!< higher order stiffness multiplier
	bool			m_btwo_pass;	//!< two-pass flag
	double			m_atol;			//!< augmentation tolerance
	double			m_gtol;			//!< gap tolerance
	double			m_stol;			//!< search tolerance
	bool			m_bsymm;		//!< use symmetric stiffness components only
	double			m_srad;			//!< contact search radius
	int				m_naugmax;		//!< maximum nr of augmentations
	int				m_naugmin;		//!< minimum nr of augmentations
	int				m_nsegup;		//!< segment update parameter
	
	double			m_epsn;			//!< normal penalty factor
	bool			m_bautopen;		//!< use autopenalty factor
	
	bool			m_btension;		//!< allow tension across interface
	
	DECLARE_PARAMETER_LIST();
};
