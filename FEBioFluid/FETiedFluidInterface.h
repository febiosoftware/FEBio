//
//  FETiedFluidInterface.hpp
//  FEBioFluid
//
//  Created by Gerard Ateshian on 8/6/18.
//  Copyright © 2018 febio.org. All rights reserved.
//

#ifndef FETiedFluidInterface_hpp
#define FETiedFluidInterface_hpp

#include "FEBioMech/FEContactInterface.h"
#include "FEBioMech/FEContactSurface.h"
#include "FEFluid.h"

//-----------------------------------------------------------------------------
class FETiedFluidSurface : public FEContactSurface
{
public:
    //! Integration point data
    class Data : public FEContactMaterialPoint
    {
    public:
        Data();
        
    public:
        vec3d   m_Gap;      //!< initial gap in reference configuration
        vec3d   m_vg;       //!< tangential velocity gap function at integration points
        vec3d   m_nu;       //!< normal at integration points
        vec2d   m_rs;       //!< natural coordinates of projection of integration point
        vec3d   m_Lmd;      //!< lagrange multipliers for tangential velocity
        vec3d   m_tv;       //!< viscous tangential traction
        double  m_Lmp;      //!< lagrange multipliers for fluid pressures
        double  m_epst;     //!< viscous traction penalty factor
        double  m_epsn;     //!< normal velocity penalty factor
        double  m_pg;       //!< pressure "gap"
        double  m_vn;       //!< normal fluid velocity gap
        FESurfaceElement*    m_pme;    //!< master element of projected integration point
    };
    
    //! constructor
    FETiedFluidSurface(FEModel* pfem);
    
    //! initialization
    bool Init();
    
    void Serialize(DumpStream& ar);
    
    //! Unpack surface element data
    void UnpackLM(FEElement& el, vector<int>& lm);

	//! create material point data
	FEMaterialPoint* CreateMaterialPoint() override;
    
public:
    void GetVelocityGap     (int nface, vec3d& vg);
    void GetPressureGap     (int nface, double& pg);
    void GetViscousTraction (int nface, vec3d& tv);
    void GetNormalVelocity  (int nface, double& vn);

protected:
    FEModel*    m_pfem;
    
public:
    int             m_dofWX, m_dofWY, m_dofWZ;
    int             m_dofEF;
};

//-----------------------------------------------------------------------------
class FETiedFluidInterface :    public FEContactInterface
{
public:
    //! constructor
    FETiedFluidInterface(FEModel* pfem);
    
    //! destructor
    ~FETiedFluidInterface();
    
    //! initialization
    bool Init() override;
    
    //! interface activation
    void Activate() override;
    
    //! serialize data to archive
    void Serialize(DumpStream& ar) override;
    
    //! return the master and slave surface
    FESurface* GetMasterSurface() override { return &m_ms; }
    FESurface* GetSlaveSurface () override { return &m_ss; }
    
    //! return integration rule class
    bool UseNodalIntegration() override { return false; }
    
    //! build the matrix profile for use in the stiffness matrix
    void BuildMatrixProfile(FEGlobalMatrix& K) override;
    
public:
    //! calculate contact forces
    void Residual(FEGlobalVector& R, const FETimeInfo& tp) override;
    
    //! calculate contact stiffness
    void StiffnessMatrix(FESolver* psolver, const FETimeInfo& tp) override;
    
    //! calculate Lagrangian augmentations
    bool Augment(int naug, const FETimeInfo& tp) override;
    
    //! update
    void Update(int niter, const FETimeInfo& tp) override;
    
protected:
    void InitialProjection(FETiedFluidSurface& ss, FETiedFluidSurface& ms);
    void ProjectSurface(FETiedFluidSurface& ss, FETiedFluidSurface& ms);
    
    //! calculate penalty factor
    void CalcAutoPressurePenalty(FETiedFluidSurface& s);
    double AutoPressurePenalty(FESurfaceElement& el, FETiedFluidSurface& s);
    
public:
    FETiedFluidSurface    m_ms;    //!< master surface
    FETiedFluidSurface    m_ss;    //!< slave surface
    
    bool            m_btwo_pass;    //!< two-pass flag
    double          m_atol;         //!< augmentation tolerance
    double          m_gtol;         //!< gap tolerance
    double          m_ptol;         //!< pressure gap tolerance
    double          m_stol;         //!< search tolerance
    double          m_srad;         //!< contact search radius
    int             m_naugmax;      //!< maximum nr of augmentations
    int             m_naugmin;      //!< minimum nr of augmentations
    
    double          m_epst;         //!< tangential viscous traction penalty factor
    double          m_epsn;         //!< normal fluid velocity penalty factor
    bool            m_bautopen;     //!< use autopenalty factor
    
    FEFluid*        m_pfluid;       //!< fluid pointer
    
    int             m_dofWX, m_dofWY, m_dofWZ;
    int             m_dofEF;
    
    DECLARE_FECORE_CLASS();
};

#endif /* FETiedFluidInterface_hpp */
