#pragma once
#include "FEContactInterface.h"
#include "FECore/FEContactSurface.h"

// Elastic sliding contact, reducing the algorithm of biphasic sliding contact
// (FESlidingInterface2) to elastic case.  The algorithm derives from Bonet
// & Wood's treatment of surface pressures

//-----------------------------------------------------------------------------
class FESlidingSurfaceBW : public FEContactSurface
{
public:
	// data for each integration point
	class Data
	{
	public:
		Data();

	public:
		double	m_gap;		//!< gap function
		double	m_Lmd;		//!< Lagrange multipliers for displacements
		double	m_Ln;		//!< net contact pressure
		double	m_epsn;		//!< penalty factor
		vec3d	m_nu;		//!< local normal
		vec2d	m_rs;		//!< natural coordinates of this integration point
		FESurfaceElement*	m_pme;	//!< projected master element
	};

public:
	//! constructor
	FESlidingSurfaceBW(FEModel* pfem);
	
	//! initialization
	bool Init();
	
	//! shallow copy
	void ShallowCopy(FESlidingSurfaceBW& s);
	
	void Serialize(DumpFile& ar);

	//! evaluate net contact force
	vec3d NetContactForce();

	//! initialize projection
	void InitProjection();
	
protected:
	FEModel*	m_pfem;
	
public:
	vector< vector<Data> >	m_Data;		//!< integration point data for all elements
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
	void ShallowCopy(FESurfacePairInteraction& ci);
	
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
