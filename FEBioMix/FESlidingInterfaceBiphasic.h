/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
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
#include "FEBioMech/FEContactInterface.h"
#include "FEBiphasicContactSurface.h"

//-----------------------------------------------------------------------------
class FEBIOMIX_API FESlidingSurfaceBiphasic : public FEBiphasicContactSurface
{
public:
    //! Integration point data
    class Data : public FEBiphasicContactPoint
    {
    public:
        Data();

		void Serialize(DumpStream& ar) override;
        
    public:
        vec3d   m_dg;       //!< vector gap
        double	m_Lmd;      //!< Lagrange multipliers for normal traction
        vec3d   m_Lmt;      //!< Lagrange multipliers for vector traction
        double	m_epsn;     //!< penalty factor
        double	m_epsp;     //!< pressure penalty factor
        double  m_p1;       //!< fluid pressure
        vec3d	m_nu;       //!< normal at integration points
        vec3d   m_s1;       //!< tangent along slip direction
        vec3d   m_tr;       //!< contact traction
        vec2d	m_rs;       //!< natural coordinates of projection
        vec2d   m_rsp;      //!< m_rs at the previous time step
        bool    m_bstick;   //!< stick flag
    };
    
public:
    //! constructor
    FESlidingSurfaceBiphasic(FEModel* pfem);
    
    //! initialization
    bool Init() override;
    
    // data serialization
    void Serialize(DumpStream& ar) override;
    
    //! initialize sliding surface and store previous values
    void InitSlidingSurface();
    
    //! evaluate net contact force
    vec3d GetContactForce() override;
    
    //! evaluate net contact area
    double GetContactArea() override;
    
    //! evaluate net fluid force
    vec3d GetFluidForce() override;
    
    //! calculate the nodal normals
    void UpdateNodeNormals();
    
    void SetPoroMode(bool bporo) { m_bporo = bporo; }
    
	//! create material point data
	FEMaterialPoint* CreateMaterialPoint() override;

public:
    void GetVectorGap      (int nface, vec3d& pg) override;
    void GetContactTraction(int nface, vec3d& pt) override;
    void GetSlipTangent    (int nface, vec3d& pt);
    void GetMuEffective    (int nface, double& pg) override;
    void GetLocalFLS       (int nface, double& pg) override;
    void GetNodalVectorGap      (int nface, vec3d* pg) override;
    void GetNodalContactPressure(int nface, double* pg) override;
    void GetNodalContactTraction(int nface, vec3d* pt) override;
    void GetStickStatus(int nface, double& pg) override;
    void EvaluateNodalContactPressures();
    void EvaluateNodalContactTractions();

private:
	void GetContactPressure(int nface, double& pg);
      
public:
    bool	m_bporo;	//!< set poro-mode
    
    vector<bool>		m_poro;	//!< surface element poro status
    vector<vec3d>		m_nn;	//!< node normals
    vector<vec3d>       m_tn;   //!< nodal contact tractions
    vector<double>      m_pn;   //!< nodal contact pressures
    
    vec3d    m_Ft;     //!< total contact force (from equivalent nodal forces)
};

//-----------------------------------------------------------------------------
class FEBIOMIX_API FESlidingInterfaceBiphasic :	public FEContactInterface
{
public:
    //! constructor
    FESlidingInterfaceBiphasic(FEModel* pfem);
    
    //! destructor
    ~FESlidingInterfaceBiphasic();
    
    //! initialization
    bool Init() override;
    
    //! interface activation
    void Activate() override;
    
    //! calculate the slip direction on the primary surface
    vec3d SlipTangent(FESlidingSurfaceBiphasic& ss, const int nel, const int nint, FESlidingSurfaceBiphasic& ms, double& dh, vec3d& r);
    
    //! calculate contact traction
    vec3d ContactTraction(FESlidingSurfaceBiphasic& ss, const int nel, const int n, FESlidingSurfaceBiphasic& ms, double& pn);
    
    //! calculate contact pressures for file output
    void UpdateContactPressures();
    
    //! serialize data to archive
    void Serialize(DumpStream& ar) override;
    
    //! mark free-draining condition
    void MarkFreeDraining();
    
    //! set free-draining condition
    void SetFreeDraining();
    
    //! return the primary and secondary surface
	FESurface* GetPrimarySurface() override { return &m_ss; }
	FESurface* GetSecondarySurface() override { return &m_ms; }

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
    void ProjectSurface(FESlidingSurfaceBiphasic& ss, FESlidingSurfaceBiphasic& ms, bool bupseg, bool bmove = false);
    
    //! calculate penalty factor
    void UpdateAutoPenalty();
    
    void CalcAutoPenalty(FESlidingSurfaceBiphasic& s);
    
    void CalcAutoPressurePenalty(FESlidingSurfaceBiphasic& s);
    double AutoPressurePenalty(FESurfaceElement& el, FESlidingSurfaceBiphasic& s);
    
public:
	FESlidingSurfaceBiphasic	m_ss;	//!< primary surface
	FESlidingSurfaceBiphasic	m_ms;	//!< secondary surface
    
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
    bool            m_bsmfls;       //!< smooth local fluid load support

    double			m_epsn;		    //!< normal penalty factor
    bool			m_bautopen;	    //!< use autopenalty factor
    bool            m_bupdtpen;     //!< update penalty at each time step
    
    double          m_mu;           //!< friction coefficient
    bool            m_bfreeze;      //!< freeze stick/slip status
    bool            m_bflips;       //!< flip primary surface normal
    bool            m_bflipm;       //!< flip secondary surface normal
    bool            m_bshellbs;     //!< flag for prescribing pressure on shell bottom for primary surface
    bool            m_bshellbm;     //!< flag for prescribing pressure on shell bottom for secondary surface

    // biphasic contact parameters
    double	        m_epsp;		    //!< flow rate penalty
    double          m_phi;          //!< solid-solid contact fraction
    
protected:
    int	m_dofP;
    
    DECLARE_FECORE_CLASS();
};
