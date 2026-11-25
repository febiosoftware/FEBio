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
#include "FEBiphasicContactSurface.h"

//-----------------------------------------------------------------------------
class FEBIOMIX_API FETiedBiphasicContactPoint: public FEBiphasicContactPoint
{
public:
    vec3d   m_Gap;  //!< initial gap (components in the current corotational frame of the primary surface)
    vec3d   m_dg;   //!< gap function at integration points
    vec3d   m_nu;   //!< normal at integration points
    vec2d   m_rs;   //!< natural coordinates of projection of integration point
    vec3d   m_Lmd;  //!< lagrange multipliers for displacements
    vec3d   m_tr;   //!< contact traction
    double  m_epsn; //!< penalty factors
    double  m_epsp; //!< pressure penalty factors

    void Init() override
    {
        FEBiphasicContactPoint::Init();
        m_Gap = m_dg = m_nu = m_Lmd = m_tr = vec3d(0,0,0);
        m_rs = vec2d(0,0);
        m_epsn = 1.0;
        m_epsp = 1.0;
    }
    
    void Serialize(DumpStream& ar) override;
};

//-----------------------------------------------------------------------------
//! The Lagrange multipliers m_Lmd are stored as spatial vectors. At the start of
//! each time step they are pushed forward with the rotation of the interface
//! (corotational frame at the primary integration point) over the previous time
//! step, so that they remain frame invariant. If binit is true, the frames are
//! only initialized.
//! The initial gap m_Gap is NOT pushed forward here: it is stored by its
//! components in the corotational frame and rotated with the frame evaluated in
//! the current configuration (see FETiedBiphasicFrameGap), since holding it fixed
//! within a time step produces spurious, step-size dependent tractions.
inline void FETiedBiphasicPushForward(FESurface& ss, bool binit)
{
    for (int i = 0; i < ss.Elements(); ++i)
    {
        FESurfaceElement& se = ss.Element(i);
        for (int j = 0; j < se.GaussPoints(); ++j)
        {
            FETiedBiphasicContactPoint& pt = static_cast<FETiedBiphasicContactPoint&>(*se.GetMaterialPoint(j));
            if (binit) pt.m_Rp = FEContactFrame();
            if (pt.m_pme == nullptr) continue;
            FEContactFrame Rn;
            Rn.Evaluate(ss, se, j, true);
            pt.m_Rp.Update(Rn, &pt.m_Lmd, 1);
        }
    }
}

//! Current image R.Gap of the initial gap stored in the corotational frame of
//! the primary surface (current configuration)
inline vec3d FETiedBiphasicFrameGap(FESurface& ss, FESurfaceElement& se, int j, const vec3d& Gap)
{
    FEContactFrame f;
    f.Evaluate(ss, se, j, false);
    return f.R(Gap);
}

//! Frame forces on the primary nodes of a tied point with an initial gap stored in
//! the corotational frame (gap function dg = (x2 - x1) - R.Gap, traction t =
//! Lm + eps*dg, nodal forces integrated over the current area, W = detJ*w).
//! These balance the moment of the traction pair (see FETiedFrameOffset).
//! fe holds 3 displacement dofs per node, primary nodes first.
inline void FETiedBiphasicFrameForce(std::vector<double>& fe, FESurface& ss, FESurfaceElement& se, int j,
    const vec3d& Gap, const vec3d& t, double eps, double W)
{
    if (Gap.norm2() == 0.0) return;
    FETiedFrameOffset fo; vec3d g1, g2; double sgn;
    fo.R.Evaluate(ss, se, j, false, g1, g2, sgn);
    fo.EvalSpatial(g1, g2, Gap, t, eps);
    double* Gr = se.Gr(j);
    double* Gs = se.Gs(j);
    for (int k = 0; k < se.Nodes(); ++k)
    {
        double T[2] = { sgn*Gr[k], sgn*Gs[k] };
        vec3d f = fo.Force(T, W);
        fe[3*k] += f.x; fe[3*k+1] += f.y; fe[3*k+2] += f.z;
    }
}

//! Tangent of the dependence of the tie on the corotational frame of the primary
//! surface (gap function dg = (x2 - x1) - R.Gap, traction t = Lm + eps*dg). This adds
//! the frame-dependent part of the Hessian of psi = Lm.dg + eps/2 dg.dg (symmetric),
//! and, if bareaterm is true, the (non-symmetric) term due to the variation of the
//! current area in the frame forces. The I-term eps*N_a*N_b and the area term of
//! the N_a*t forces are evaluated by the caller.
//! ndpn is the number of dofs per node in ke (displacements first).
template <class MAT>
inline void FETiedBiphasicFrameStiffness(MAT& ke, FESurface& ss, FESurfaceElement& se, int j,
    const vec3d& Gap, const vec3d& t, double eps, double w, const double* Hs, int nseln, const double* Hm, int nmeln,
    int ndpn, bool bareaterm)
{
    if (Gap.norm2() == 0.0) return;
    FETiedFrameOffset fo; vec3d g1, g2; double sgn;
    fo.R.Evaluate(ss, se, j, false, g1, g2, sgn);
    fo.EvalSpatial(g1, g2, Gap, t, eps);
    vec3d gc[2];
    ss.CoBaseVectors(se, j, gc);
    double W = (gc[0] ^ gc[1]).norm()*w;
    double* Gr = se.Gr(j);
    double* Gs = se.Gs(j);
    const int MN = FEElement::MAX_NODES;
    int off[2*MN]; double Y[2*MN]; double T[2*MN][2];
    for (int a = 0; a < nseln; ++a) { off[a] = ndpn*a; Y[a] = -Hs[a]; T[a][0] = sgn*Gr[a]; T[a][1] = sgn*Gs[a]; }
    for (int b = 0; b < nmeln; ++b) { off[nseln+b] = ndpn*(nseln+b); Y[nseln+b] = Hm[b]; T[nseln+b][0] = T[nseln+b][1] = 0.0; }
    FEContactFrame::AddHessian(ke, nseln + nmeln, off, Y, T, eps, false, fo.Hyg, fo.Hgg, W);
    if (bareaterm) fo.AddAreaTangent(ke, nseln, off, T, gc[0], gc[1], Gr, Gs, w, 1.0);
}


//-----------------------------------------------------------------------------
class FEBIOMIX_API FETiedBiphasicSurface : public FEBiphasicContactSurface
{
public:
	//! constructor
	FETiedBiphasicSurface(FEModel* pfem);
	
	//! initialization
	bool Init() override;
	
	//! create material point data
	FEMaterialPoint* CreateMaterialPoint() override;

	//! calculate the nodal normals
	void UpdateNodeNormals();
	
	void Serialize(DumpStream& ar) override;
	
	void SetPoroMode(bool bporo) { m_bporo = bporo; }
	
public:
    void GetVectorGap      (int nface, vec3d& pg) override;
    void GetContactTraction(int nface, vec3d& pt) override;
    
public:
	bool				m_bporo;	//!< set poro-mode
	
    vector<bool>			m_poro;	//!< surface element poro status
    vector<vec3d>			m_nn;	//!< node normals
};

//-----------------------------------------------------------------------------
class FEBIOMIX_API FETiedBiphasicInterface :	public FEContactInterface
{
public:
	//! constructor
	FETiedBiphasicInterface(FEModel* pfem);
	
	//! destructor
	~FETiedBiphasicInterface();
	
	//! initialization
	bool Init() override;
	
	//! interface activation
	void Activate() override;

	//! push stored vectors forward with the interface rotation (frame invariance)
	void PrepStep() override;

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

protected:
	void InitialProjection(FETiedBiphasicSurface& ss, FETiedBiphasicSurface& ms);
	void ProjectSurface(FETiedBiphasicSurface& ss, FETiedBiphasicSurface& ms);
	
	//! calculate penalty factor
    void UpdateAutoPenalty();
    
	void CalcAutoPenalty(FETiedBiphasicSurface& s);
	
	void CalcAutoPressurePenalty(FETiedBiphasicSurface& s);
	double AutoPressurePenalty(FESurfaceElement& el, FETiedBiphasicSurface& s);
	
public:
	FETiedBiphasicSurface	m_ss;	//!< primary surface
	FETiedBiphasicSurface	m_ms;	//!< secondary surface
	
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
	
	double			m_epsn;			//!< normal penalty factor
	bool			m_bautopen;		//!< use autopenalty factor
    bool            m_bupdtpen;     //!< update penalty at each time step

	// biphasic contact parameters
	double			m_epsp;		//!< flow rate penalty

protected:
	int	m_dofP;
	
	DECLARE_FECORE_CLASS();
};
