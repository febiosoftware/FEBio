#pragma once
#include "FEContactInterface.h"
#include "FEContactSurface.h"

// Elastic sliding contact, reducing the algorithm of biphasic sliding contact
// (FESlidingInterface2) to elastic case.  The algorithm derives from Bonet
// & Wood's treatment of surface pressures
// Modified by Brandon Zimmerman & Gerard Ateshian on 5/28/17 to include friction

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
        double	m_gap;		//!< normal gap function
        vec3d   m_dg;       //!< vector gap
        double	m_Lmd;		//!< Lagrange multiplier for normal traction
        double	m_Ln;		//!< net contact pressure
        double	m_epsn;		//!< penalty factor
        vec3d   m_Lmt;      //!< Lagrange multipliers for vector traction
        vec3d	m_nu;		//!< local normal
        vec3d   m_s1;       //!< tangent along slip direction
        vec3d   m_tr;       //!< contact traction
        vec2d	m_rs;		//!< natural coordinates of this integration point
        vec2d   m_rsp;      //!< m_rs at the previous time step
        bool    m_bstick;   //!< stick flag
        FESurfaceElement*	m_pme;	//!< projected master element
        FESurfaceElement*   m_pmep; //!< m_pme at the previous time step
    };
    
public:
    //! constructor
    FESlidingSurfaceBW(FEModel* pfem);
    
    //! initialization
    bool Init();
    
    void Serialize(DumpStream& ar);
    
    //! initialize sliding surface and store previous values
    void InitSlidingSurface();
    
    //! evaluate net contact force
    vec3d GetContactForce();
    
    //! evaluate net contact area
    double GetContactArea();
    
public:
    void GetContactGap     (int nface, double& pg);
    void GetVectorGap      (int nface, vec3d& pg);
    void GetContactPressure(int nface, double& pg);
    void GetContactTraction(int nface, vec3d& pt);
    void GetNodalContactGap     (int nface, double* pg);
    void GetNodalVectorGap      (int nface, vec3d* pg);
    void GetNodalContactPressure(int nface, double* pg);
    void GetNodalContactTraction(int nface, vec3d* pt);
    void GetStickStatus(int nface, double& pg);
    
protected:
    FEModel*	m_pfem;
    
public:
    vector< vector<Data> >	m_Data;		//!< integration point data for all elements
    
    vec3d    m_Ft;     //!< total contact force (from equivalent nodal forces)
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
    
    //! calculate the slip direction on the primary surface
    vec3d SlipTangent(FESlidingSurfaceBW& ss, const int nel, const int nint, FESlidingSurfaceBW& ms, double& dh, vec3d& r);
    
    //! calculate contact traction
    vec3d ContactTraction(FESlidingSurfaceBW& ss, const int nel, const int n, FESlidingSurfaceBW& ms, double& pn);
   
    //! calculate contact pressures for file output
    void UpdateContactPressures();
    
    //! serialize data to archive
    void Serialize(DumpStream& ar);
    
    //! return the master and slave surface
    FESurface* GetMasterSurface() { return &m_ms; }
    FESurface* GetSlaveSurface () { return &m_ss; }
    
    //! return integration rule class
    bool UseNodalIntegration() { return false; }
    
    //! build the matrix profile for use in the stiffness matrix
    void BuildMatrixProfile(FEGlobalMatrix& K);

public:
	//! calculate contact forces
	void Residual(FEGlobalVector& R, const FETimeInfo& tp);

	//! calculate contact stiffness
	void StiffnessMatrix(FESolver* psolver, const FETimeInfo& tp);

	//! calculate Lagrangian augmentations
	bool Augment(int naug, const FETimeInfo& tp);

	//! update
	void Update(int niter, const FETimeInfo& tp);

protected:
    void ProjectSurface(FESlidingSurfaceBW& ss, FESlidingSurfaceBW& ms, bool bupseg, bool bmove = false);
    
    //! calculate penalty factor
    void CalcAutoPenalty(FESlidingSurfaceBW& s);
    
public:
    FESlidingSurfaceBW	m_ms;	//!< master surface
    FESlidingSurfaceBW	m_ss;	//!< slave surface
    
    int				m_knmult;		//!< higher order stiffness multiplier
    bool			m_btwo_pass;	//!< two-pass flag
    double			m_atol;			//!< augmentation tolerance
    double			m_gtol;			//!< normal gap tolerance
    double			m_stol;			//!< search tolerance
    bool			m_bsymm;		//!< use symmetric stiffness components only
    double			m_srad;			//!< contact search radius
    int				m_naugmax;		//!< maximum nr of augmentations
    int				m_naugmin;		//!< minimum nr of augmentations
    int				m_nsegup;		//!< segment update parameter
    bool			m_breloc;		//!< node relocation on activation
    bool            m_bsmaug;       //!< smooth augmentation
    
    double			m_epsn;			//!< normal penalty factor
    bool			m_bautopen;		//!< use autopenalty factor
    
    bool			m_btension;		//!< allow tension across interface
    
    double          m_mu;           //!< friction coefficient
    
    bool            m_bfreeze;      //!< freeze stick/slip status
    
    DECLARE_PARAMETER_LIST();
};
