//
//  FETiedElasticInterface.hpp
//  FEBioMech
//
//  Created by Gerard Ateshian on 2/28/18.
//  Copyright © 2018 febio.org. All rights reserved.
//

#ifndef FETiedElasticInterface_hpp
#define FETiedElasticInterface_hpp

#include "FEBioMech/FEContactInterface.h"
#include "FEContactSurface.h"

//-----------------------------------------------------------------------------
class FETiedElasticSurface : public FEContactSurface
{
public:
    //! Integration point data
    class Data
    {
    public:
        Data();
        
    public:
        vec3d    m_Gap;     //!< initial gap in reference configuration
        vec3d    m_dg;      //!< gap function at integration points
        vec3d    m_nu;      //!< normal at integration points
        vec2d    m_rs;      //!< natural coordinates of projection of integration point
        vec3d    m_Lmd;     //!< lagrange multipliers for displacements
        vec3d    m_tr;      //!< contact traction
        double   m_epsn;    //!< penalty factors
        FESurfaceElement*    m_pme;    //!< master element of projected integration point
    };
    
    //! constructor
    FETiedElasticSurface(FEModel* pfem);
    
    //! initialization
    bool Init();
    
    //! calculate the nodal normals
    void UpdateNodeNormals();
    
    void Serialize(DumpStream& ar);
    
public:
    void GetVectorGap      (int nface, vec3d& pg);
    void GetContactTraction(int nface, vec3d& pt);
    
protected:
    FEModel*    m_pfem;
    
public:
    vector< vector<Data> >  m_Data; //!< integration point data
    vector<vec3d>           m_nn;   //!< node normals
};

//-----------------------------------------------------------------------------
class FETiedElasticInterface : public FEContactInterface
{
public:
    //! constructor
    FETiedElasticInterface(FEModel* pfem);
    
    //! destructor
    ~FETiedElasticInterface();
    
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
    void InitialProjection(FETiedElasticSurface& ss, FETiedElasticSurface& ms);
    void ProjectSurface(FETiedElasticSurface& ss, FETiedElasticSurface& ms);
    
    //! calculate penalty factor
    void CalcAutoPenalty(FETiedElasticSurface& s);
    
public:
    FETiedElasticSurface    m_ms;    //!< master surface
    FETiedElasticSurface    m_ss;    //!< slave surface
    
    int         m_knmult;       //!< higher order stiffness multiplier
    bool        m_btwo_pass;    //!< two-pass flag
    double      m_atol;         //!< augmentation tolerance
    double      m_gtol;         //!< gap tolerance
    double      m_stol;         //!< search tolerance
    bool        m_bsymm;        //!< use symmetric stiffness components only
    double      m_srad;         //!< contact search radius
    int         m_naugmax;      //!< maximum nr of augmentations
    int         m_naugmin;      //!< minimum nr of augmentations
    
    double      m_epsn;         //!< normal penalty factor
    bool        m_bautopen;     //!< use autopenalty factor

    DECLARE_FECORE_CLASS();
};

#endif /* FETiedElasticInterface_hpp */
