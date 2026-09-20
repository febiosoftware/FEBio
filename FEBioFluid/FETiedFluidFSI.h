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
//! Tied fluid-FSI interface, implementing the formulation of Section 7.9 of the
//! FEBio Theory Manual.
//!
//! This interface follows the tied fluid interface of Section 7.8, with two
//! modifications. First, the solid displacement u is added to the list of nodal
//! degrees of freedom, so each node carries seven dofs (u, w, J) instead of four.
//! Second, the nodal fluid degree of freedom in the FSI framework is not the fluid
//! velocity v^f but the fluid velocity relative to the solid, w = v^f - v^s.
//! Because the interface also enforces [[u]] = 0, and hence [[v^s]] = 0,
//! constraining [[w]] = 0 is equivalent to the no-slip condition [[v^f]] = 0.
//!
//! Three constraints are enforced at each integration point, each with its own
//! penalty factor (eq. (7.9.4)-(7.9.6)):
//!
//!     t^s   = lambda_s + eps_s*g     solid traction    (vectorial gap g)
//!     t^tau = lambda_t + eps_t*g^w   viscous traction  (relative velocity gap g^w)
//!     w_n   = lambda_p + eps_n*pi    normal velocity   (dilatation gap pi)
//!
//! All three gaps are oriented as (secondary - primary), i.e. g^w = w(2) - w(1)
//! of eq. (7.9.4) and pi = J(2) - J(1). This orientation is the one for which the
//! stiffness blocks of eqs. (7.9.13)-(7.9.16) hold as written.
//!
//! ---------------------------------------------------------------------------------
//! FSI-SOLID PAIRS (an extension beyond Section 7.9)
//!
//! Section 7.9 ties two FSI domains: every term of eq. (7.9.3) pairs dv^s, dw and dJ
//! across both surfaces. This class additionally supports tying an FSI domain to a
//! plain elastic solid domain, whose nodes carry only u. That case is NOT symmetric
//! and needs a different constraint set:
//!
//!   1. [[u]] = 0, enforced exactly as in Section 7.9 with the eps_s penalty. This
//!      already transmits action and reaction between the two displacement fields.
//!   2. w.n = 0 on the FSI side: the wall is impermeable, but tangential slip is left
//!      free. This is one-sided, and is obtained by treating the solid wall as
//!      w(2) = 0 and then keeping only the normal part of the resulting gap:
//!
//!          g^w = -(w.n) n ,   t^tau = (lambda.n + eps_t (g^w.n)) n
//!
//!      so the traction is purely normal. Note this is weaker than the full no-slip
//!      condition w = 0 that the User Manual prescribes for a conforming FSI-solid
//!      boundary: the fluid may slide along the wall and develops no boundary layer
//!      there. The solid still feels the viscous shear, because the fluid traction
//!      transferred in item 4 below is the full sigma^f.n.
//!
//!      Because t^tau now depends on the normal, and hence on u, its linearization
//!      carries a dn/du term that the Section 7.9 blocks do not have. It is included.
//!      It sits in the w-row / u-column block, whose transpose counterpart is
//!      identically zero, so like the K^fs and k^fs terms of eq. (7.9.14) it cannot be
//!      symmetrized and symmetric_stiffness drops it. (The solid A_c term is different:
//!      it has a diagonal counterpart, so that one is symmetrized rather than dropped.)
//!      Dropping it costs the quadratic convergence rate but not the converged
//!      solution, since none of these terms touches the residual.
//!   3. No dilatation constraint: J is free at an FSI-solid boundary, as it is in the
//!      conforming case. The pi gap and its blocks are dropped.
//!   4. The fluid traction t^f = sigma^f.n must be transferred to the solid. The
//!      penalty tie on its own only balances t^s(1) against t^solid(2); t^f appears in
//!      neither u-equation, which is exactly why a conforming mesh needs a separate
//!      "fluid-FSI traction" load. This class therefore owns an FEFluidFSITraction on
//!      the FSI surface and adds its contribution, so the tie carries t^f across.
//!      In the matched-shape-function limit this reproduces the conforming force
//!      balance F^s_int(1) + F^solid_int(2) = T^f.
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
                            //!< surface, -(w.n) n against a solid wall
        vec3d    m_wr;      //!< relative fluid velocity w on the primary surface. Needed
                            //!< by the tangent of the normal-only wall constraint, whose
                            //!< traction depends on u through the normal.
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
    //! The alpha_f prefactors of eqs. (7.9.13)-(7.9.16) come entirely from
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
    //! (returns nullptr if the element is not backed by a fluid-FSI material)
    FEFluid* GetFluid(FESurfaceElement& el);

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
        TIE_FSI_FSI   = 0,  //!< two FSI domains, the formulation of Section 7.9
        TIE_FSI_SOLID = 1   //!< FSI primary against a plain solid secondary
    };

public:
    FETiedFluidFSISurface    m_ss;    //!< primary surface
    FETiedFluidFSISurface    m_ms;    //!< secondary surface

    int         m_knmult;       //!< higher order stiffness multiplier
    bool        m_btwo_pass;    //!< two-pass flag
    double      m_atol;         //!< augmentation tolerance
    double      m_gtol;         //!< solid displacement gap tolerance
    double      m_wtol;         //!< relative fluid velocity gap tolerance. Against a
                                //!< solid wall this bounds |w.n| only: the tangential
                                //!< slip is unconstrained and is not measured.
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
    //! dofs of the FSI surface, from where the displacement tie carries it to the
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
