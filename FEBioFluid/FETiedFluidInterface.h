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
#include <FEBioMech/FEContactInterface.h>
#include <FEBioMech/FEContactSurface.h>
#include <FECore/FESolidElement.h>
#include "FEFluidMaterial.h"

//-----------------------------------------------------------------------------
class FEBIOFLUID_API FETiedFluidSurface : public FEContactSurface
{
public:
    //! Integration point data
    class Data : public FEContactMaterialPoint
    {
    public:
        Data();

		void Serialize(DumpStream& ar) override;
        
    public:
        vec3d   m_Gap;      //!< initial gap in reference configuration
        vec3d   m_vg;       //!< velocity gap function at integration points
        vec3d   m_nu;       //!< normal at integration points
        vec2d   m_rs;       //!< natural coordinates of projection of integration point
        vec3d   m_Lmd;      //!< lagrange multipliers for velocity
        vec3d   m_tv;       //!< viscous traction
        double  m_Lmp;      //!< lagrange multipliers for fluid dilatations
        double  m_epst;     //!< viscous traction penalty factor
        double  m_epsn;     //!< normal velocity penalty factor
        double  m_Jg;       //!< dilatation "gap", pi = J(2) - J(1)
        double  m_vn;       //!< normal velocity
    };
    
    //! constructor
    FETiedFluidSurface(FEModel* pfem);
    
    //! initialization
    bool Init() override;
    
    //! Unpack surface element data
    void UnpackLM(FEElement& el, vector<int>& lm) override;

	//! create material point data
	FEMaterialPoint* CreateMaterialPoint() override;
    
public:
    void GetVelocityGap     (int nface, vec3d& vg);
    void GetDilatationGap   (int nface, double& Jg);
    void GetViscousTraction (int nface, vec3d& tv);
    void GetNormalVelocity  (int nface, double& vn);
    double GetArea          (FESurfaceElement& el, bool breference = false);
    double GetVolume        (FESolidElement& el);
   
public:
	FEDofList	m_dofWE;
};

//-----------------------------------------------------------------------------
class FEBIOFLUID_API FETiedFluidInterface :    public FEContactInterface
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
    
    //! return the primary and secondary surfaces
	FESurface* GetPrimarySurface() override { return &m_s1; }
	FESurface* GetSecondarySurface() override { return &m_s2; }
    
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
    //! initial projection; bfirst is true only on the first (primary->secondary) pass,
    //! so that the m_bfreedofs option only ever frees dofs on the secondary surface
    void InitialProjection(FETiedFluidSurface& ss, FETiedFluidSurface& ms, bool bfirst);
    void ProjectSurface(FETiedFluidSurface& ss, FETiedFluidSurface& ms);

    //! Record a dof that the m_bfreedofs option has to keep open, and open it.
    //! nodeIndex is an index into the mesh node array, not a node ID.
    void RecordFreeDof(int nodeIndex, int dof);

    //! Re-open every dof recorded by RecordFreeDof. Run at the start of every analysis
    //! step, because a step activates its own BCs and those may re-constrain these dofs.
    void ReleaseFreeDofs();

    //! Hook ReleaseFreeDofs() up to CB_STEP_ACTIVE. Idempotent, and called both from
    //! Init() and from Serialize() on load, since a dump restart does not re-run Init().
    void RegisterFreeDofsCallback();

    //! CB_STEP_ACTIVE callback: re-applies ReleaseFreeDofs() after the step's boundary
    //! conditions have been activated and before the solver numbers the equations.
    static bool free_dofs_cb(FEModel* pfem, unsigned int nwhen, void* pd);
    
    //! return the fluid material shared by all elements attached to this surface
    //! (returns nullptr if the surface is not backed by a single fluid material)
    FEFluidMaterial* GetFluidMaterial(FETiedFluidSurface& s);
    
    //! calculate penalty factor
    void CalcAutoViscousTractionPenalty(FETiedFluidSurface& s);
    void CalcAutoNormalVelocityPenalty(FETiedFluidSurface& s);
    
    double AutoViscousTractionPenalty(FESurfaceElement& el, FETiedFluidSurface& s);
    double AutoNormalVelocityPenalty(FESurfaceElement& el, FETiedFluidSurface& s);
public:
	FETiedFluidSurface    m_s1;    //!< primary surface
	FETiedFluidSurface    m_s2;    //!< secondary surface
    
    bool            m_btwo_pass;    //!< two-pass flag
    double          m_atol;         //!< augmentation tolerance
    double          m_gtol;         //!< velocity gap tolerance
    double          m_etol;         //!< dilatation gap tolerance
    double          m_stol;         //!< search tolerance
    double          m_srad;         //!< contact search radius
    int             m_naugmax;      //!< maximum nr of augmentations
    int             m_naugmin;      //!< minimum nr of augmentations
    
    double          m_epst;         //!< tangential viscous traction penalty factor
    double          m_epsn;         //!< normal fluid velocity penalty factor
    bool            m_bautopen;     //!< use autopenalty factor
    
    bool            m_bfreedofs;    //!< flag to free constrained/fixed DOFS on secondary surface

    //! dofs released by the m_bfreedofs option, as parallel (node index, dof) lists.
    //! Kept so the release can be re-applied at every step boundary.
    std::vector<int>    m_freeNode;
    std::vector<int>    m_freeDof;
    bool                m_bfreecb = false;  //!< true once the CB_STEP_ACTIVE callback is registered

    FEFluidMaterial* m_pfluid = nullptr;    //!< fluid pointer (set in Init)

	FEDofList		m_dofWE;
   
    DECLARE_FECORE_CLASS();
};
