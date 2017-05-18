//
//  FESlidingInterfaceBiphasic.hpp
//  FEBioMix
//
//  Created by Gerard Ateshian on 5/1/16.
//  Copyright © 2016 febio.org. All rights reserved.
//

#ifndef FESlidingInterfaceBiphasic_hpp
#define FESlidingInterfaceBiphasic_hpp

#include "FEBioMech/FEContactInterface.h"
#include "FEBiphasicContactSurface.h"

//-----------------------------------------------------------------------------
class FESlidingSurfaceBiphasic : public FEBiphasicContactSurface
{
public:
    //! Integration point data
    class Data
    {
    public:
        Data();
        
    public:
        double	m_gap;      //!< gap function
        vec3d   m_dg;       //!< vector gap
        double	m_Lmd;      //!< Lagrange multipliers for normal traction
        vec3d   m_Lmt;      //!< Lagrange multipliers for vector traction
        double	m_Lmp;      //!< lagrange multipliers for fluid pressures
        double	m_Ln;       //!< net contact pressure
        double	m_epsn;     //!< penalty factor
        double	m_epsp;     //!< pressure penatly factor
        double	m_pg;       //!< pressure "gap" for biphasic contact
        double  m_p1;       //!< fluid pressure
        vec3d	m_nu;       //!< normal at integration points
        vec3d   m_s1;       //!< tangent along slip direction
        vec3d   m_tr;       //!< contact traction
        vec2d	m_rs;       //!< natrual coordinates of projection
        vec2d   m_rsp;      //!< m_rs at the previous time step
        bool    m_bstick;   //!< stick flag
        FESurfaceElement*	m_pme;	//!< master element
        FESurfaceElement*   m_pmep; //!< m_pme at the previous time step
    };
    
public:
    //! constructor
    FESlidingSurfaceBiphasic(FEModel* pfem);
    
    //! initialization
    bool Init();
    
    // data serialization
    void Serialize(DumpStream& ar);
    
    //! initialize sliding surface and store previous values
    void InitSlidingSurface();
    
    //! evaluate net contact force
    vec3d GetContactForce();
    
    //! evaluate net contact area
    double GetContactArea();
    
    //! evaluate net fluid force
    vec3d GetFluidForce();
    
    //! calculate the nodal normals
    void UpdateNodeNormals();
    
    void SetPoroMode(bool bporo) { m_bporo = bporo; }
    
public:
    void GetContactGap     (int nface, double& pg);
    void GetContactPressure(int nface, double& pg);
    void GetContactTraction(int nface, vec3d& pt);
    void GetNodalContactGap     (int nface, double* pg);
    void GetNodalContactPressure(int nface, double* pg);
    void GetNodalContactTraction(int nface, vec3d* pt);
    void GetNodalPressureGap    (int nface, double* pg);
    void EvaluateNodalContactTractions();
    
protected:
    FEModel*	m_pfem;
    
public:
    bool	m_bporo;	//!< set poro-mode
    
    vector< vector<Data> >	m_Data;	//!< integration point data
    vector<bool>		m_poro;	//!< surface element poro status
    vector<vec3d>		m_nn;	//!< node normals
    vector<vec3d>       m_tn;   //!< nodal contact tractions

    vec3d    m_Ft;     //!< total contact force (from equivalent nodal forces)
};

//-----------------------------------------------------------------------------
class FESlidingInterfaceBiphasic :	public FEContactInterface
{
public:
    //! constructor
    FESlidingInterfaceBiphasic(FEModel* pfem);
    
    //! destructor
    ~FESlidingInterfaceBiphasic();
    
    //! initialization
    bool Init();
    
    //! interface activation
    void Activate();
    
    //! update
    void Update(int niter);
    
    //! calculate the slip direction on the primary surface
    vec3d SlipTangent(FESlidingSurfaceBiphasic& ss, const int nel, const int nint, FESlidingSurfaceBiphasic& ms, double& dh, vec3d& r);
    
    //! calculate contact traction
    vec3d ContactTraction(FESlidingSurfaceBiphasic& ss, const int nel, const int n, FESlidingSurfaceBiphasic& ms, double& pn);
    
    //! calculate contact forces
    void ContactForces(FEGlobalVector& R);
    
    //! calculate contact stiffness
    void ContactStiffness(FESolver* psolver);
    
    //! calculate contact pressures for file output
    void UpdateContactPressures();
    
    //! calculate Lagrangian augmentations
    bool Augment(int naug);
    
    //! serialize data to archive
    void Serialize(DumpStream& ar);
    
    //! mark free-draining condition
    void MarkFreeDraining();
    
    //! set free-draining condition
    void SetFreeDraining();
    
    //! return the master and slave surface
    FESurface* GetMasterSurface() { return &m_ms; }
    FESurface* GetSlaveSurface () { return &m_ss; }
    
    //! return integration rule class
    bool UseNodalIntegration() { return false; }
    
    //! build the matrix profile for use in the stiffness matrix
    void BuildMatrixProfile(FEGlobalMatrix& K);
    
protected:
    void ProjectSurface(FESlidingSurfaceBiphasic& ss, FESlidingSurfaceBiphasic& ms, bool bupseg, bool bmove = false);
    
    //! calculate penalty factor
    void CalcAutoPenalty(FESlidingSurfaceBiphasic& s);
    
    void CalcAutoPressurePenalty(FESlidingSurfaceBiphasic& s);
    double AutoPressurePenalty(FESurfaceElement& el, FESlidingSurfaceBiphasic& s);
    
public:
    FESlidingSurfaceBiphasic	m_ms;	//!< master surface
    FESlidingSurfaceBiphasic	m_ss;	//!< slave surface
    
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
    bool			m_breloc;		//!< node relocation on startup
    bool            m_bsmaug;       //!< smooth augmentation
    
    double			m_epsn;		//!< normal penalty factor
    bool			m_bautopen;	//!< use autopenalty factor
    
    // biphasic contact parameters
    double	m_epsp;		//!< flow rate penalty
    
    double          m_mu;           //!< friction coefficient
    double          m_phi;          //!< solid-solid contact fraction
    
protected:
    int	m_dofP;
    
    DECLARE_PARAMETER_LIST();
};

#endif /* FESlidingInterfaceBiphasic_hpp */
