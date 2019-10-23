/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, The Trustees of Columbia University in 
the City of New York, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/



#pragma once
#include "FEContactInterface.h"
#include "FEContactSurface.h"

// Elastic sliding contact, reducing the algorithm of biphasic sliding contact
// (FESlidingInterface2) to elastic case.  The algorithm derives from Bonet
// & Wood's treatment of surface pressures

//-----------------------------------------------------------------------------
class FESlidingSurfaceBW : public FEContactSurface
{
public:
    // data for each integration point
    class Data : public FEContactMaterialPoint
    {
    public:
        Data();
        
    public:
        vec3d   m_dg;       //!< vector gap
        double	m_Lmd;		//!< Lagrange multiplier for normal traction
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
    bool Init() override;
    
    void Serialize(DumpStream& ar) override;
    
    //! initialize sliding surface and store previous values
    void InitSlidingSurface();
    
    //! evaluate net contact force
    vec3d GetContactForce() override;
    
    //! evaluate net contact area
    double GetContactArea() override;
    
	//! create material point data
	FEMaterialPoint* CreateMaterialPoint() override;

public:
    void GetVectorGap      (int nface, vec3d& pg) override;
    void GetContactTraction(int nface, vec3d& pt) override;
    void GetNodalVectorGap      (int nface, vec3d* pg) override;
    void GetNodalContactPressure(int nface, double* pg) override;
    void GetNodalContactTraction(int nface, vec3d* pt) override;
    void GetStickStatus(int nface, double& pg) override;
    
protected:
    FEModel*	m_pfem;
    
public:
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
    bool Init() override;
    
    //! interface activation
    void Activate() override;
    
    //! calculate the slip direction on the primary surface
    vec3d SlipTangent(FESlidingSurfaceBW& ss, const int nel, const int nint, FESlidingSurfaceBW& ms, double& dh, vec3d& r);
    
    //! calculate contact traction
    vec3d ContactTraction(FESlidingSurfaceBW& ss, const int nel, const int n, FESlidingSurfaceBW& ms, double& pn);
   
    //! calculate contact pressures for file output
    void UpdateContactPressures();
    
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
	void LoadVector(FEGlobalVector& R, const FETimeInfo& tp) override;

	//! calculate contact stiffness
	void StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp) override;

	//! calculate Lagrangian augmentations
	bool Augment(int naug, const FETimeInfo& tp) override;

	//! update
	void Update() override;

protected:
    void ProjectSurface(FESlidingSurfaceBW& ss, FESlidingSurfaceBW& ms, bool bupseg, bool bmove = false);
    
    //! calculate penalty factor
    void CalcAutoPenalty(FESlidingSurfaceBW& s);

	void SerializePointers(FESlidingSurfaceBW& ss, FESlidingSurfaceBW& ms, DumpStream& ar);
    
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
    bool            m_bflipm;       //!< flip master normal
    bool            m_bflips;       //!< flip slave normal

    DECLARE_FECORE_CLASS();
};
