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
#include "FEFluid.h"
#include "FEFluidFSI.h"
#include "FEFluidFSITraction.h"

//-----------------------------------------------------------------------------
//! Tied fluid-FSI interface, implementing the formulation of Section 7.11 of the
//! FEBio Theory Manual.
//!
//! This interface extends the tied fluid interface of Section 7.9, with two
//! modifications. First, the solid displacement u is added to the list of nodal
//! degrees of freedom, so each node carries seven dofs (u, w, J) instead of four.
//! Second, the nodal fluid degree of freedom in the FSI framework is not the fluid
//! velocity v^f but the fluid velocity relative to the solid, w = v^f - v^s.
//! Because the interface also enforces [[u]] = 0, and hence [[v^s]] = 0,
//! constraining [[w]] = 0 is equivalent to the no-slip condition [[v^f]] = 0.
//!
//! Three constraints are enforced at each integration point, each with its own
//! penalty factor (eqs. (7.11.7)-(7.11.8), augmented as in eq. (7.6.5)):
//!
//!     t^s   = lambda_s + eps_s*g     solid traction    (vectorial gap g)
//!     t^tau = lambda_t + eps_t*g^w   viscous traction  (relative velocity gap g^w)
//!     w_n   = lambda_p + eps_n*pi    normal velocity   (dilatation gap pi)
//!
//! All three gaps are oriented as (secondary - primary), eqs. (7.11.4)-(7.11.6):
//! g = x(2) - x(1), g^w = w(2) - w(1) and pi = J(2) - J(1). This orientation is the
//! one for which the stiffness blocks of eqs. (7.11.23)-(7.11.25) hold as written.
//!
//! ---------------------------------------------------------------------------------
//! FSI-SOLID PAIRS (Section 7.10, Tied Fluid-FSI-Solid Traction Interface)
//!
//! Section 7.11 ties two FSI domains: every term of eq. (7.11.1) pairs dv^s, dw and
//! dJ across both surfaces. This class additionally supports tying an FSI domain to
//! a plain elastic solid domain, whose nodes carry only u. Following Section 7.10,
//! the constraint set in that case is:
//!
//!   1. [[u]] = 0, the tied-elastic contact of Section 7.6 with the eps_s penalty.
//!      This transmits action and reaction between the two displacement fields.
//!   2. The no-slip condition w = 0 on the FSI surface Gamma(1), which Section 7.10
//!      requires to be prescribed in addition to the tie. It is enforced here by the
//!      viscous traction penalty, treating the solid wall as w(2) = 0:
//!
//!          g^w = -w(1) ,   t^tau = lambda_t + eps_t g^w
//!
//!      so the w-w block is the isotropic block of eq. (7.11.23) (K^ww_ac), and the
//!      only u-dependence of the w-rows is through J_eta (K^wu_ac).
//!   3. No dilatation constraint: J is free at an FSI-solid boundary, as it is in the
//!      conforming case. The pi gap and its blocks are dropped.
//!   4. The fluid traction t^f = sigma^f.n must be transferred to the solid. This is
//!      the term dF = -int dv^s(1).t^f(1) da(1) of Section 7.10, which is identical
//!      to the fluid-FSI traction of Section 3.6.11 applied on Gamma(1) with its
//!      original sign (n(1) outward from the FSI domain). The tie on its own only
//!      balances t^s(1) against t^s(2); this class therefore owns an
//!      FEFluidFSITraction on the FSI surface and adds its contribution, so that
//!      t^s(2) ~ -t^f(1) is transmitted to the solid.
//!
//! In this mode the FSI surface must be the primary surface, two_pass is not allowed
//! (the constraint set is not symmetric), and only FEFluidFSI is supported: biphasic-
//! and multiphasic-FSI need the solid-volume-fraction weighting of
//! FEBiphasicFSITraction and are rejected rather than silently mistreated.

//-----------------------------------------------------------------------------
class FEBIOFLUID_API FETiedFluidFSISurface : public FEContactSurface
{
public:
    //! Integration point data
    class Data : public FEContactMaterialPoint
    {
    public:
        Data();

		void Serialize(DumpStream& ar) override;

        void Init() override;

    public:
        vec3d    m_Gap;     //!< initial gap in reference configuration
        vec3d    m_dg;      //!< vectorial gap g at integration points
        vec3d    m_gw;      //!< relative fluid velocity gap: w(2) - w(1) against an FSI
                            //!< surface, -w(1) against a solid wall (no-slip, Sec. 7.10)
        double   m_Jg;      //!< dilatation gap pi = J(2) - J(1)
        vec3d    m_nu;      //!< normal at integration points
        vec2d    m_rs;      //!< natural coordinates of projection of integration point
        vec3d    m_Lmd;     //!< Lagrange multiplier lambda_s for the solid traction
        vec3d    m_Lmt;     //!< Lagrange multiplier lambda_t for the viscous traction
        double   m_Lmp;     //!< Lagrange multiplier lambda_p for the normal velocity
        vec3d    m_ts;      //!< solid traction t^s
        vec3d    m_tv;      //!< viscous traction t^tau
        double   m_wn;      //!< normal relative fluid velocity w_n
        double   m_epss;    //!< solid traction penalty factor
        double   m_epst;    //!< viscous traction penalty factor
        double   m_epsn;    //!< dilatation penalty factor
    };

    //! constructor
    FETiedFluidFSISurface(FEModel* pfem);

    //! initialization
    bool Init() override;

    //! calculate the nodal normals
    void UpdateNodeNormals();

    void Serialize(DumpStream& ar) override;

	//! create material point data
	FEMaterialPoint* CreateMaterialPoint() override;

    //! Unpack the LM vector. Seven dofs per node: u (3), w (3) and J (1).
    void UnpackLM(FEElement& el, vector<int>& lm) override;

    //! surface area of a face, and volume of the solid element backing it.
    //! Their ratio is the element thickness used by the automatic penalties.
    double GetArea          (FESurfaceElement& el, bool breference = false);
    double GetVolume        (FESolidElement& el);

    //! Geometry in the alpha_f-weighted intermediate configuration.
    //!
    //! The alpha_f prefactors of eqs. (7.11.9)-(7.11.15) come entirely from
    //! D(x_alpha)[Du] = alpha_f Du, so those blocks are the consistent tangent only if
    //! the position, the normal and the covariant basis vectors (hence J_eta) are all
    //! evaluated at the intermediate configuration rather than at the end of the step.
    //! This is also the convention used throughout the rest of FEBioFluid's FSI loads.
    vec3d Local2GlobalAlpha (FESurfaceElement& el, int n, double alpha);
    vec3d Local2GlobalAlpha (FESurfaceElement& el, double r, double s, double alpha);
    void  CoBaseVectorsAlpha(FESurfaceElement& el, int n, double alpha, vec3d g[2]);
    vec3d SurfaceNormalAlpha(FESurfaceElement& el, int n, double alpha);

public:
    void GetVectorGap      (int nface, vec3d& pg) override;
    void GetContactTraction(int nface, vec3d& pt) override;

    //! evaluate net contact force
    vec3d GetContactForce() override;

    //! evaluate net contact area
    double GetContactArea() override;

    //! Does this surface belong to an FSI domain (u, w, J) or a plain solid (u only)?
    //! Set by the interface once it has identified the backing materials.
    void SetFSI(bool b) { m_bfsi = b; }
    bool IsFSI() const { return m_bfsi; }

    //! Dofs per node in the LM vector. This is ALWAYS 7, on a solid surface too.
    //!
    //! A plain solid node carries no w and no J, but FEResidualVector and
    //! FERigidSolver::RigidStiffness recover the nodal stride as ndof/en.size() and so
    //! assume one uniform stride across the whole element vector. Packing 7 for the
    //! FSI face and 3 for the solid face would make that quotient a truncated average
    //! and index past the end of the node list. Instead a solid surface emits -1 for
    //! its w and J entries, which every assembler already skips, exactly as it would
    //! for any other inactive dof (FESolver::InitEquations leaves them at -1).
    int DofsPerNode() const { return 7; }

public:
    vector<vec3d>           m_nn;   //!< node normals
    vector<vec3d>           m_Fn;   //!< nodal forces

    //! the seven dofs of this surface, in the order u (3), w (3), J (1).
    //! A solid surface resolves all seven but only ever unpacks the first three.
    FEDofList               m_dofUWE;

protected:
    bool                    m_bfsi = true;  //!< true if backed by an FSI material
};

//-----------------------------------------------------------------------------
class FEBIOFLUID_API FETiedFluidFSI : public FEContactInterface
{
public:
    //! constructor
    FETiedFluidFSI(FEModel* pfem);

    //! destructor
    ~FETiedFluidFSI();

    //! initialization
    bool Init() override;

    //! interface activation
    void Activate() override;

    //! serialize data to archive
    void Serialize(DumpStream& ar) override;

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

    //! called at the start of each time step; recomputes the automatic penalties
    //! when update_penalty is set
    void PrepStep() override;

protected:
    //! initial projection. bfirst is true only on the first (primary->secondary)
    //! pass, so that one-sided bookkeeping is never repeated on the swapped pass.
    void InitialProjection(FETiedFluidFSISurface& ss, FETiedFluidFSISurface& ms, bool bfirst);
    void ProjectSurface(FETiedFluidFSISurface& ss, FETiedFluidFSISurface& ms);

    //! recalculate all automatic penalty factors
    void UpdateAutoPenalty();

    //! solid traction penalty eps_s (elastic auto penalty of the backing solid)
    void CalcAutoPenalty(FETiedFluidFSISurface& s);

    //! viscous traction penalty eps_t (fluid viscosity over element thickness)
    void CalcAutoViscousTractionPenalty(FETiedFluidFSISurface& s);
    double AutoViscousTractionPenalty(FESurfaceElement& el, FETiedFluidFSISurface& s);

    //! dilatation penalty eps_n (viscous relaxation scaling, else acoustic impedance)
    void CalcAutoNormalVelocityPenalty(FETiedFluidFSISurface& s);
    double AutoNormalVelocityPenalty(FESurfaceElement& el, FETiedFluidFSISurface& s);

    //! return the fluid of the FSI material backing this surface element
    //! (returns nullptr if the element is not backed by a fluid-FSI material).
    //! If ppe is given, it receives the FSI solid element the fluid was taken from,
    //! which is the element whose material points and volume must be used.
    FEFluid* GetFluid(FESurfaceElement& el, FESolidElement** ppe = nullptr);

    //! Identify the materials backing each surface and set m_mode accordingly.
    //! Returns false, with a logged error, on any unsupported combination.
    bool DetectMode();

    //! Is this surface backed by an FSI material? Classification is per face, by the
    //! element the face belongs to. bbiphasic is set if a biphasic- or multiphasic-FSI
    //! material was found, which this class rejects; bmixed is set if the surface
    //! straddles FSI and non-FSI domains, which cannot be classified at all.
    bool IsFSISurface(FETiedFluidFSISurface& s, bool& bbiphasic, bool& bmixed);

public:
    //! which kind of pair this interface ties
    enum TieMode {
        TIE_FSI_FSI   = 0,  //!< two FSI domains, the formulation of Section 7.11
        TIE_FSI_SOLID = 1   //!< FSI primary against a plain solid secondary (Section 7.10)
    };

public:
    FETiedFluidFSISurface    m_ss;    //!< primary surface
    FETiedFluidFSISurface    m_ms;    //!< secondary surface

    int         m_knmult;       //!< higher order stiffness multiplier
    bool        m_btwo_pass;    //!< two-pass flag
    double      m_atol;         //!< augmentation tolerance
    double      m_gtol;         //!< solid displacement gap tolerance
    double      m_wtol;         //!< relative fluid velocity gap tolerance. Against a
                                //!< solid wall this bounds |w| (no-slip, Sec. 7.10).
    double      m_etol;         //!< dilatation gap tolerance
    double      m_stol;         //!< search tolerance
    bool        m_bsymm;        //!< use symmetric stiffness components only
    double      m_srad;         //!< contact search radius
    int         m_naugmax;      //!< maximum nr of augmentations
    int         m_naugmin;      //!< minimum nr of augmentations

    double      m_epss;         //!< solid traction penalty factor
    double      m_epst;         //!< viscous traction penalty factor
    double      m_epsn;         //!< dilatation penalty factor
    bool        m_bautopen;     //!< use autopenalty factor
    bool        m_bupdtpen;     //!< update penalty at each time step

	bool            m_bflips;       //!< flip primary surface normal
	bool            m_bflipm;       //!< flip secondary surface normal

protected:
    int             m_mode = TIE_FSI_FSI;  //!< set by DetectMode

    //! In TIE_FSI_SOLID mode this transfers the fluid traction sigma^f.n onto the u
    //! dofs of the FSI surface (the dF term of Section 7.10), from where the displacement tie carries it to the
    //! solid. Owned by this interface; null in TIE_FSI_FSI mode, where both sides
    //! are FSI and no transfer is needed.
    FEFluidFSITraction*  m_pFSItrac = nullptr;

protected:
    // degrees of freedom
    FEDofList   m_dofU, m_dofSU, m_dofW;
    int         m_dofEF;

protected:
    bool                m_bshellb;  //!< flag for prescribing traction on shell bottom

    DECLARE_FECORE_CLASS();
};
