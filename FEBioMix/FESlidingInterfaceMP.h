#pragma once
#include "FEBioMech/FEContactInterface.h"
#include "FECore/FEContactSurface.h"

//-----------------------------------------------------------------------------
class FESlidingSurfaceMP : public FEContactSurface
{
public:
	//! integration point data
	class Data
	{
	public:
		Data();
        
	public:
		double	m_gap;	//!< gap function at integration points
		double	m_Lmd;	//!< Lagrange multipliers for displacements
		double	m_Lmp;  //!< Lagrange multipliers for fluid pressure
		double	m_Ln;	//!< net contact pressure
		double	m_epsn;	//!< displacement penalty factors
		double	m_epsp;	//!< pressure penalty factors
		double	m_pg;	//!< pressure "gap"
		vec3d	m_nu;	//!< normal at integration points
		vec2d	m_rs;	//!< natural coordinates of projection of integration point
		vector<double>	m_Lmc;	//!< Lagrange multipliers for solute concentrations
		vector<double>	m_epsc;	//!< concentration penatly factors
		vector<double>	m_cg;	//!< concentration "gap"
		FESurfaceElement*	m_pme;	//!< master element of projected integration point
	};
    
public:
	//! constructor
	FESlidingSurfaceMP(FEModel* pfem);
	
	//! destructor
	~FESlidingSurfaceMP() {}
	
	//! initialization
	bool Init();
	
	//! shallow copy
	void ShallowCopy(DumpStream& dmp, bool bsave);
	
	//! evaluate net contact force
	vec3d GetContactForce();
	
	//! evaluate net fluid force
	vec3d GetFluidForce();
	
	//! calculate the nodal normals
	void UpdateNodeNormals();
	
	void Serialize(DumpFile& ar);
	
	void SetPoroMode(bool bporo) { m_bporo = bporo; }
	
public:
	void GetNodalContactGap     (int nface, double* pg);
	void GetNodalContactPressure(int nface, double* pg);
	void GetNodalContactTraction(int nface, vec3d* tn);
	
protected:
	FEModel*	m_pfem;
	
public:
	vector< vector<Data> >		m_Data; //!< integration point data
    
	bool						m_bporo;	//!< set poro-mode
	bool						m_bsolu;	//!< set solute-mode
	
	vector<vec3d>				m_nn;	//!< node normals
	
	vector<int>					m_sid;	//!< list of solute id's for this surface
};

//-----------------------------------------------------------------------------
class FESlidingInterfaceMP : public FEContactInterface
{
public:
	//! constructor
	FESlidingInterfaceMP(FEModel* pfem);
	
	//! destructor
	~FESlidingInterfaceMP();
	
	//! initialization
	bool Init();
	
	//! interface activation
	void Activate();

	//! update
	void Update(int niter);
	
	//! Create a shallow copy
	void ShallowCopy(DumpStream& dmp, bool bsave);
	
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
	
	//! mark ambient condition 
	void MarkAmbient();
	
	//! set ambient condition 
	void SetAmbient();
	
	//! return the master and slave surface
	FESurface* GetMasterSurface() { return &m_ms; }
	FESurface* GetSlaveSurface () { return &m_ss; }
    
	//! return integration rule class
	bool UseNodalIntegration() { return false; }
    
	//! build the matrix profile for use in the stiffness matrix
	void BuildMatrixProfile(FEStiffnessMatrix& K);
    
protected:
	void ProjectSurface(FESlidingSurfaceMP& ss, FESlidingSurfaceMP& ms, bool bupseg);
	
	//! calculate penalty factor
	void CalcAutoPenalty(FESlidingSurfaceMP& s);
	
	void CalcAutoPressurePenalty(FESlidingSurfaceMP& s);
	double AutoPressurePenalty(FESurfaceElement& el, FESlidingSurfaceMP& s);
	
	void CalcAutoConcentrationPenalty(FESlidingSurfaceMP& s, const int isol);
	double AutoConcentrationPenalty(FESurfaceElement& el, FESlidingSurfaceMP& s, const int isol);
	
public:
	FESlidingSurfaceMP	m_ms;	//!< master surface
	FESlidingSurfaceMP	m_ss;	//!< slave surface
	
	int				m_knmult;		//!< higher order stiffness multiplier
	bool			m_btwo_pass;	//!< two-pass flag
	double			m_atol;			//!< augmentation tolerance
	double			m_gtol;			//!< gap tolerance
	double			m_ptol;			//!< pressure gap tolerance
	double			m_ctol;			//!< concentration gap tolerance
	double			m_stol;			//!< search tolerance
	bool			m_bsymm;		//!< use symmetric stiffness components only
	double			m_srad;			//!< contact search radius
	int				m_naugmax;		//!< maximum nr of augmentations
	int				m_naugmin;		//!< minimum nr of augmentations
	int				m_nsegup;		//!< segment update parameter
	
	double			m_epsn;		//!< normal penalty factor
	bool			m_bautopen;	//!< use autopenalty factor
	
	// multiphasic contact parameters
	double	m_epsp;					//!< fluid volumetric flow rate penalty
	double	m_epsc;					//!< solute molar flow rate penalty
	double	m_Rgas;					//!< universal gas constant
	double	m_Tabs;					//!< absolute temperature
	double	m_ambp;					//!< ambient pressure
	double	m_ambc[MAX_CDOFS];		//!< ambient concentration
	vector<int> m_sid;				//!< list of solute ids common to both contact surfaces
	vector<int> m_ssl;				//!< list of slave surface solutes common to both contact surfaces
	vector<int> m_msl;				//!< list of master surface solutes common to both contact surfaces
	
	DECLARE_PARAMETER_LIST();
};
