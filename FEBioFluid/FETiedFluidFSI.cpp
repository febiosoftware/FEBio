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



#include "stdafx.h"
#include "FETiedFluidFSI.h"
#include "FEFluid.h"
#include "FEFluidFSI.h"
#include "FEBiphasicFSI.h"
#include "FEFluidFSITraction.h"
#include "FEBioFSI.h"
#include <FECore/FEModel.h>
#include <FECore/FEAnalysis.h>
#include <FECore/FENormalProjection.h>
#include <FECore/FELinearSystem.h>
#include <FECore/FESolidElement.h>
#include <FECore/log.h>

//-----------------------------------------------------------------------------
// Define tied fluid-FSI interface parameters
BEGIN_FECORE_CLASS(FETiedFluidFSI, FEContactInterface)
	ADD_PARAMETER(m_laugon   , "laugon"             )->setLongName("Enforcement method")->setEnums("PENALTY\0AUGLAG\0");
	ADD_PARAMETER(m_atol     , "tolerance"          );
	ADD_PARAMETER(m_gtol     , "gaptol"             );
	ADD_PARAMETER(m_wtol     , "vtol"               );
	ADD_PARAMETER(m_etol     , "etol"               );
	ADD_PARAMETER(m_epss     , "solid_penalty"      );
	ADD_PARAMETER(m_epst     , "traction_penalty"   );
	ADD_PARAMETER(m_epsn     , "dilatation_penalty" );
	ADD_PARAMETER(m_bautopen , "auto_penalty"       );
    ADD_PARAMETER(m_bupdtpen , "update_penalty"     );
	ADD_PARAMETER(m_btwo_pass, "two_pass"           );
	ADD_PARAMETER(m_knmult   , "knmult"             );
	ADD_PARAMETER(m_stol     , "search_tol"         );
	ADD_PARAMETER(m_bsymm    , "symmetric_stiffness");
	ADD_PARAMETER(m_srad     , "search_radius"      );
	ADD_PARAMETER(m_naugmin  , "minaug"             );
	ADD_PARAMETER(m_naugmax  , "maxaug"             );
	ADD_PARAMETER(m_bflips   , "flip_primary"       );
	ADD_PARAMETER(m_bflipm   , "flip_secondary"     );
    ADD_PARAMETER(m_bshellb  , "shell_bottom");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FETiedFluidFSISurface::Data::Data()
{
    m_Gap = vec3d(0,0,0);
    m_dg  = vec3d(0,0,0);
    m_gw  = vec3d(0,0,0);
    m_wr  = vec3d(0,0,0);
    m_Jg  = 0.0;
    m_nu  = vec3d(0,0,0);
    m_rs  = vec2d(0,0);
    m_Lmd = vec3d(0,0,0);
    m_Lmt = vec3d(0,0,0);
    m_Lmp = 0.0;
    m_ts  = vec3d(0,0,0);
    m_tv  = vec3d(0,0,0);
    m_wn  = 0.0;
    m_epss= 1.0;
    m_epst= 1.0;
    m_epsn= 1.0;
}

//-----------------------------------------------------------------------------
void FETiedFluidFSISurface::Data::Serialize(DumpStream& ar)
{
	FEContactMaterialPoint::Serialize(ar);
	ar & m_Gap;
	ar & m_dg;
	ar & m_gw;
	ar & m_wr;
	ar & m_Jg;
	ar & m_nu;
	ar & m_rs;
	ar & m_Lmd;
	ar & m_Lmt;
	ar & m_Lmp;
	ar & m_ts;
	ar & m_tv;
	ar & m_wn;
	ar & m_epss;
	ar & m_epst;
	ar & m_epsn;
}

void FETiedFluidFSISurface::Data::Init()
{
    FEContactMaterialPoint::Init();
    m_Gap = vec3d(0, 0, 0);
    m_dg  = vec3d(0, 0, 0);
    m_gw  = vec3d(0, 0, 0);
    m_wr  = vec3d(0, 0, 0);
    m_Jg  = 0.0;
    m_nu  = vec3d(0, 0, 0);
    m_rs  = vec2d(0, 0);
    m_Lmd = vec3d(0, 0, 0);
    m_Lmt = vec3d(0, 0, 0);
    m_Lmp = 0.0;
    m_ts  = vec3d(0, 0, 0);
    m_tv  = vec3d(0, 0, 0);
    m_wn  = 0.0;
    m_epss= 1.0;
    m_epst= 1.0;
    m_epsn= 1.0;
}

//-----------------------------------------------------------------------------
// FETiedFluidFSISurface
//-----------------------------------------------------------------------------

FETiedFluidFSISurface::FETiedFluidFSISurface(FEModel* pfem) : FEContactSurface(pfem), m_dofUWE(pfem)
{
}

//-----------------------------------------------------------------------------
bool FETiedFluidFSISurface::Init()
{
    // initialize surface data first
    if (FEContactSurface::Init() == false) return false;

    // Set the dof list: seven dofs per node, in the order u (3), w (3), J (1).
    // The order in which the variables are added here is the dof ordering used by
    // UnpackLM, and therefore the ordering of the element vectors and matrices.
    //
    // NOTE: this guard matters. FESurfaceLoad::Init calls Init() on its surface, so
    // when the interface hands this surface to the fluid traction load it owns, we
    // land here a second time. FEDofList::AddVariable appends unconditionally, so
    // without the guard the list would grow to 14 and the size check below would fail.
    if (m_dofUWE.Size() == 0)
    {
        if (m_dofUWE.AddVariable(FEBioFSI::GetVariableName(FEBioFSI::DISPLACEMENT)) == false) return false;
        if (m_dofUWE.AddVariable(FEBioFSI::GetVariableName(FEBioFSI::RELATIVE_FLUID_VELOCITY)) == false) return false;
        if (m_dofUWE.AddVariable(FEBioFSI::GetVariableName(FEBioFSI::FLUID_DILATATION)) == false) return false;
    }
    if (m_dofUWE.Size() != 7) return false;

    // allocate node normals
    m_nn.assign(Nodes(), vec3d(0,0,0));

    // initialize nodal force vector
    m_Fn.assign(Nodes(), vec3d(0, 0, 0));

    return true;
}

//-----------------------------------------------------------------------------
//! Unpack the LM vector: seven entries per node, ordered u, w, J.
//!
//! On a plain solid surface the w and J entries are emitted as -1, which every
//! assembler skips. Keeping the stride uniform at 7 is deliberate: see DofsPerNode.
void FETiedFluidFSISurface::UnpackLM(FEElement& el, vector<int>& lm)
{
    int N = el.Nodes();
    const int ndpn = 7;
    const int nact = (m_bfsi ? 7 : 3);
    lm.resize(N*ndpn);
    for (int i=0; i<N; ++i)
    {
        int n = el.m_node[i];
        FENode& node = m_pMesh->Node(n);
        vector<int>& id = node.m_ID;

        for (int k=0; k<nact; ++k) lm[ndpn*i+k] = id[m_dofUWE[k]];
        for (int k=nact; k<ndpn; ++k) lm[ndpn*i+k] = -1;
    }
}

//-----------------------------------------------------------------------------
//! This function calculates the node normal. Due to the piecewise continuity
//! of the surface elements this normal is not uniquely defined so in order to
//! obtain a unique normal the normal is averaged for each node over all the
//! element normals at the node

void FETiedFluidFSISurface::UpdateNodeNormals()
{
    int N = Nodes(), i, j, ne, jp1, jm1;
    vec3d y[FEElement::MAX_NODES], n;

    // zero nodal normals
    zero(m_nn);

    // loop over all elements
    for (i=0; i<Elements(); ++i)
    {
        FESurfaceElement& el = Element(i);
        ne = el.Nodes();

        // get the nodal coordinates
        for (j=0; j<ne; ++j) y[j] = Node(el.m_lnode[j]).m_rt;

        // calculate the normals
        for (j=0; j<ne; ++j)
        {
            jp1 = (j+1)%ne;
            jm1 = (j+ne-1)%ne;
            n = (y[jp1] - y[j]) ^ (y[jm1] - y[j]);
            m_nn[el.m_lnode[j]] += n;
        }
    }

    // normalize all vectors
    for (i=0; i<N; ++i) m_nn[i].unit();
}

//-----------------------------------------------------------------------------
//! create material point data
FEMaterialPoint* FETiedFluidFSISurface::CreateMaterialPoint()
{
	return new FETiedFluidFSISurface::Data;
}

//-----------------------------------------------------------------------------
void FETiedFluidFSISurface::Serialize(DumpStream& ar)
{
	FEContactSurface::Serialize(ar);
	ar & m_nn & m_Fn;
    if (ar.IsShallow()) return;
    ar & m_dofUWE;
    // m_bfsi decides how many dofs this surface unpacks and, on the secondary side,
    // whether the dilatation is constrained at all. A dump restart does not re-run
    // Init, so losing it would silently change the physics after the restart.
    ar & m_bfsi;
}

//-----------------------------------------------------------------------------
//! Surface area of a face, in the reference or the current configuration.
double FETiedFluidFSISurface::GetArea(FESurfaceElement& el, bool breference)
{
    int ni = el.GaussPoints();
    vec3d rt[FEElement::MAX_NODES];
    if (breference) {
        for (int i=0; i<el.Nodes(); ++i) rt[i] = m_pMesh->Node(el.m_node[i]).m_r0;
    }
    else
        GetNodalCoordinates(el, 1.0, rt);
    double* gw = el.GaussWeights();
    double area = 0;
    for (int i=0; i<ni;++i) {
        FEMaterialPoint& mp = *el.GetMaterialPoint(i);
        vec3d dxr = el.eval_deriv1(rt, mp.m_index);
        vec3d dxs = el.eval_deriv2(rt, mp.m_index);
        // normal and area element
        vec3d n = dxr ^ dxs;
        double da = n.unit();
        area += da*gw[i];
    }
    return area;
}

//-----------------------------------------------------------------------------
//! Reference volume of a solid element.
double FETiedFluidFSISurface::GetVolume(FESolidElement& el)
{
    int ni = el.GaussPoints();
    double* gw = el.GaussWeights();
    double volume = 0;
    for (int i=0; i<ni; ++i) {
        double detJ = 1./(el.m_J0i[i].det());
        volume += detJ*gw[i];
    }
    return volume;
}

//-----------------------------------------------------------------------------
//! position of integration point n in the intermediate configuration
vec3d FETiedFluidFSISurface::Local2GlobalAlpha(FESurfaceElement& el, int n, double alpha)
{
    vec3d rt[FEElement::MAX_NODES];
    GetNodalCoordinates(el, alpha, rt);
    return el.eval(rt, n);
}

//-----------------------------------------------------------------------------
//! position at parametric coordinates (r,s) in the intermediate configuration
vec3d FETiedFluidFSISurface::Local2GlobalAlpha(FESurfaceElement& el, double r, double s, double alpha)
{
    vec3d rt[FEElement::MAX_NODES];
    GetNodalCoordinates(el, alpha, rt);
    return el.eval(rt, r, s);
}

//-----------------------------------------------------------------------------
//! covariant basis vectors g_1, g_2 at integration point n, intermediate configuration
void FETiedFluidFSISurface::CoBaseVectorsAlpha(FESurfaceElement& el, int n, double alpha, vec3d g[2])
{
    vec3d rt[FEElement::MAX_NODES];
    GetNodalCoordinates(el, alpha, rt);
    g[0] = el.eval_deriv1(rt, n);
    g[1] = el.eval_deriv2(rt, n);
}

//-----------------------------------------------------------------------------
//! unit normal at integration point n, intermediate configuration
vec3d FETiedFluidFSISurface::SurfaceNormalAlpha(FESurfaceElement& el, int n, double alpha)
{
    vec3d g[2];
    CoBaseVectorsAlpha(el, n, alpha, g);
    vec3d nu = g[0] ^ g[1];
    nu.unit();
    return nu;
}

//-----------------------------------------------------------------------------
void FETiedFluidFSISurface::GetVectorGap(int nface, vec3d& pg)
{
    FESurfaceElement& el = Element(nface);
    int ni = el.GaussPoints();
    pg = vec3d(0,0,0);
	for (int k = 0; k < ni; ++k)
	{
		Data& data = static_cast<Data&>(*el.GetMaterialPoint(k));
		pg += data.m_dg;
	}
    pg /= ni;
}

//-----------------------------------------------------------------------------
//! The reported contact traction is the solid traction t^s, which is the part of
//! the mixture traction of eq. (7.9.1) that this interface transmits through the
//! solid skeleton.
void FETiedFluidFSISurface::GetContactTraction(int nface, vec3d& pt)
{
    FESurfaceElement& el = Element(nface);
    int ni = el.GaussPoints();
    pt = vec3d(0,0,0);
	for (int k = 0; k < ni; ++k)
	{
		Data& data = static_cast<Data&>(*el.GetMaterialPoint(k));
		pt += data.m_ts;
	}
    pt /= ni;
}

//-----------------------------------------------------------------------------
vec3d FETiedFluidFSISurface::GetContactForce()
{
    // initialize contact force
    vec3d f(0, 0, 0);

    // loop over all elements of the primary surface
    for (size_t i = 0; i < m_Fn.size(); ++i) f += m_Fn[i];

    return f;
}

//-----------------------------------------------------------------------------
//! evaluate net contact area
double FETiedFluidFSISurface::GetContactArea()
{
    // initialize contact area
    double a = 0;

    // loop over all elements of the primary surface
    for (int n = 0; n < Elements(); ++n)
    {
        FESurfaceElement& el = Element(n);
        int nint = el.GaussPoints();

        // evaluate the contact force for that element
        for (int i = 0; i < nint; ++i)
        {
            // get data for this integration point
            Data& data = static_cast<Data&>(*el.GetMaterialPoint(i));
            double T = data.m_ts.norm2();
            if (data.m_pme && (T != 0.0))
            {
                // get the base vectors
                vec3d g[2];
                CoBaseVectors(el, i, g);

                // normal (magnitude = area)
                vec3d n = g[0] ^ g[1];

                // gauss weight
                double w = el.GaussWeights()[i];

                // contact force
                a += n.norm() * w;
            }
        }
    }

    return a;
}

//-----------------------------------------------------------------------------
// FETiedFluidFSI
//-----------------------------------------------------------------------------

FETiedFluidFSI::FETiedFluidFSI(FEModel* pfem) : FEContactInterface(pfem), m_ss(pfem), m_ms(pfem), m_dofU(pfem), m_dofSU(pfem), m_dofW(pfem)
{
    static int count = 1;
    SetID(count++);

    // initial values
    m_knmult = 1;
    m_atol = 0.1;
    m_epss = 1;
    m_epst = 1;
    m_epsn = 1;
    m_stol = 0.01;
    m_bsymm = false;   // the geometric blocks of eq. (7.9.14) are not symmetric
    m_srad = 1.0;
    m_gtol = -1;    // we use augmentation tolerance by default
    m_wtol = -1;    // we use augmentation tolerance by default
    m_etol = -1;    // we use augmentation tolerance by default
    m_bautopen = false;
    m_bupdtpen = false;
    m_btwo_pass = false;

    m_naugmin = 0;
    m_naugmax = 10;

	m_bflips = false;
	m_bflipm = false;

    // set parents
    m_ss.SetContactInterface(this);
    m_ms.SetContactInterface(this);

    m_ss.SetSibling(&m_ms);
    m_ms.SetSibling(&m_ss);

    m_bshellb = false;

    m_dofEF = -1;

    // get the degrees of freedom
    if (pfem)
    {
        m_dofU.AddVariable(FEBioFSI::GetVariableName(FEBioFSI::DISPLACEMENT));
        m_dofSU.AddVariable(FEBioFSI::GetVariableName(FEBioFSI::SHELL_DISPLACEMENT));
        m_dofW.AddVariable(FEBioFSI::GetVariableName(FEBioFSI::RELATIVE_FLUID_VELOCITY));
        m_dofEF = GetDOFIndex(FEBioFSI::GetVariableName(FEBioFSI::FLUID_DILATATION), 0);
    }
}

//-----------------------------------------------------------------------------

FETiedFluidFSI::~FETiedFluidFSI()
{
    delete m_pFSItrac;
    m_pFSItrac = nullptr;
}

//-----------------------------------------------------------------------------
//! Is this surface backed by an FSI material? Sets bbiphasic if a biphasic- or
//! multiphasic-FSI material is found; those derive from FEFluidFSI but need the
//! solid-volume-fraction weighting of FEBiphasicFSITraction, which this class does
//! not implement, so they are reported rather than silently treated as fluid-FSI.
bool FETiedFluidFSI::IsFSISurface(FETiedFluidFSISurface& s, bool& bbiphasic, bool& bmixed)
{
    bbiphasic = false;
    bmixed = false;
    int nfsi = 0, nsolid = 0;

    for (int i=0; i<s.Elements(); ++i)
    {
        FESurfaceElement& el = s.Element(i);

        // Classify by the element this face actually belongs to, NOT by "either
        // neighbour". At a conforming FSI|solid boundary FindElements populates both
        // m_elem[0] and m_elem[1] on BOTH surfaces, so an either-neighbour test would
        // report the solid surface as FSI and silently select the FSI-FSI mode: no
        // fluid traction transfer, plus a spurious dilatation constraint reading J
        // from solid nodes that do not have it. The owning element is the one the
        // face normal points away from, which is what FacePointing returns +1 for.
        FEElement* own = el.m_elem[0].pe;
        for (int k=0; k<2; ++k)
        {
            FEElement* pe = el.m_elem[k].pe;
            if (pe == nullptr) continue;
            if (s.FacePointing(el, *pe) > 0) { own = pe; break; }
        }
        if (own == nullptr) { ++nsolid; continue; }

        FEMaterial* pmat = GetFEModel()->GetMaterial(own->GetMatID());
        // check the biphasic variants first: FEBiphasicFSI derives from FEFluidFSI
        if (dynamic_cast<FEBiphasicFSI*>(pmat)) { bbiphasic = true; ++nfsi; }
        else if (dynamic_cast<FEFluidFSI*>(pmat)) ++nfsi;
        else ++nsolid;
    }

    // A surface backed by a mix of FSI and solid elements cannot be classified, and
    // treating it as either would silently drop the constraints on the other faces.
    // There is no correct answer here, so this is an error rather than a warning.
    if ((nfsi > 0) && (nsolid > 0)) {
        bmixed = true;
        feLogError("Tied fluid-FSI interface %d: surface \"%s\" has %d face(s) backed by an FSI\n"
                   "material and %d backed by something else. Each surface of this interface must be\n"
                   "entirely one or the other, otherwise the constraints on the remaining faces would\n"
                   "be silently dropped. Split the surface.", GetID(), s.GetName().c_str(), nfsi, nsolid);
    }

    return (nfsi > 0);
}

//-----------------------------------------------------------------------------
//! Identify what kind of pair this interface ties, and reject the combinations
//! this class does not support.
bool FETiedFluidFSI::DetectMode()
{
    bool bbp1 = false, bbp2 = false, bmx1 = false, bmx2 = false;
    bool bfsi1 = IsFSISurface(m_ss, bbp1, bmx1);
    bool bfsi2 = IsFSISurface(m_ms, bbp2, bmx2);

    if (bmx1 || bmx2) return false;

    if (bbp1 || bbp2) {
        feLogError("Tied fluid-FSI interface %d: biphasic-FSI and multiphasic-FSI materials are\n"
                   "not supported. Their interface traction carries a solid volume fraction\n"
                   "weighting (see FEBiphasicFSITraction) that this interface does not apply.", GetID());
        return false;
    }

    if (bfsi1 && bfsi2)
    {
        // the formulation of Section 7.9
        m_mode = TIE_FSI_FSI;
        m_ss.SetFSI(true);
        m_ms.SetFSI(true);
        feLog(" tied fluid-FSI interface # %d: FSI-FSI pair (Section 7.9)\n", GetID());
        return true;
    }

    if (bfsi1 && !bfsi2)
    {
        m_mode = TIE_FSI_SOLID;
        m_ss.SetFSI(true);
        m_ms.SetFSI(false);

        if (m_btwo_pass) {
            feLogError("Tied fluid-FSI interface %d: two_pass cannot be used with an FSI-solid pair.\n"
                       "The constraint set is not symmetric: w.n = 0 and the fluid traction transfer\n"
                       "apply to the FSI surface only, so the two surfaces cannot exchange roles.", GetID());
            return false;
        }

        feLog(" tied fluid-FSI interface # %d: FSI-solid pair\n", GetID());
        feLog("    the fluid traction is transferred to the solid by this interface;\n");
        feLog("    do not also define a fluid-FSI traction on the primary surface.\n");
        return true;
    }

    if (!bfsi1 && bfsi2) {
        feLogError("Tied fluid-FSI interface %d: the FSI surface must be the PRIMARY surface.\n"
                   "This interface found a solid primary and an FSI secondary; swap them in the\n"
                   "contact definition. (The impermeability constraint and the fluid traction are both\n"
                   "integrated over the FSI surface, which is the surface this interface loops over.)", GetID());
        return false;
    }

    feLogError("Tied fluid-FSI interface %d: neither surface is backed by a fluid-FSI material.\n"
               "This interface requires at least one FSI domain.", GetID());
    return false;
}

//-----------------------------------------------------------------------------
bool FETiedFluidFSI::Init()
{
    // Init can run more than once (parameter optimization, FEModel::Reset). Drop any
    // load left from a previous pass before the mode is redetermined, so a run that
    // now resolves to an FSI-FSI pair cannot keep contributing a stale traction.
    delete m_pFSItrac;
    m_pFSItrac = nullptr;

    // The dof lists are resolved in the constructor, which cannot report an error.
    // Check them here: without the FSI variables this interface has nothing to tie,
    // and ProjectSurface would index past the end of an empty dof list.
    if ((m_dofU.Size() != 3) || (m_dofW.Size() != 3) || (m_dofEF < 0)) {
        feLogError("Tied fluid-FSI interface %d requires the solid displacement, relative fluid\n"
                   "velocity and fluid dilatation degrees of freedom. This interface can only be\n"
                   "used in a fluid-FSI analysis.", GetID());
        return false;
    }

    // Ask the surfaces to locate both of their adjacent elements. The FSI surface of an
    // FSI-solid pair is typically an exterior face with a single neighbour, which is the
    // case FEFluidFSITraction::Activate handles through its bself branch, but the flag
    // also has to be set before FESurface::Init for that lookup to happen at all.
    m_ss.SetInterfaceStatus(true);
    m_ms.SetInterfaceStatus(true);

    // initialize surface data
    if (m_ss.Init() == false) return false;
    if (m_ms.Init() == false) return false;

    // Identify what we are tying. This has to follow the surface Init above, which is
    // what resolves el.m_elem[], and precede anything that depends on the dof counts.
    if (DetectMode() == false) return false;
    // NOTE: shell_bottom is applied to the primary surface only. With two_pass the two
    // surfaces exchange roles, so the flag would have to apply to whichever surface is
    // currently primary; rather than guess, refuse the combination.
    m_ss.SetShellBottom(m_bshellb);
    if (m_bshellb && m_btwo_pass) {
        feLogError("Tied fluid-FSI interface %d: shell_bottom cannot be combined with\n"
                   "two_pass, because the flag applies to the primary surface only.", GetID());
        return false;
    }

	// Flip secondary and primary surfaces, if requested.
	// Note that we turn off those flags because otherwise we keep flipping, each time we get here (e.g. in optimization)
	if (m_bflips) { m_ss.Invert(); m_bflips = false; }
	if (m_bflipm) { m_ms.Invert(); m_bflipm = false; }

    // In an FSI-solid pair the fluid traction has to be transferred to the solid, and
    // the penalty tie alone does not do it: t^f appears in neither u-equation. Own an
    // FEFluidFSITraction on the FSI surface and add its contribution, which is what a
    // conforming mesh achieves with a separately declared "fluid-FSI traction" load.
    if (m_mode == TIE_FSI_SOLID)
    {
        FEModel* fem = GetFEModel();
        m_pFSItrac = new FEFluidFSITraction(fem);
        m_pFSItrac->SetSurface(&m_ss);
        // The flag has to be set on the LOAD, not just on the surface: it also selects
        // whether the traction is assembled onto the displacement or the shell
        // displacement dofs. Setting it here also makes the load's Init apply it to
        // the surface, so it must precede that call.
        m_pFSItrac->SetShellBottom(m_bshellb);
        if (m_pFSItrac->Init() == false) {
            feLogError("Tied fluid-FSI interface %d: could not initialize the fluid traction\n"
                       "transfer on the primary surface.", GetID());
            return false;
        }
    }

    return true;
}

//-----------------------------------------------------------------------------
//! build the matrix profile for use in the stiffness matrix
//!
//! Each node of a tied pair couples all seven of its dofs to all seven dofs of
//! every node of the opposing face: the solid blocks of eq. (7.9.13) couple u to
//! u, the blocks of eq. (7.9.14) couple w and J to u, and eqs. (7.9.15)-(7.9.16)
//! couple w to w and J to J. A single lm vector holding all seven dofs of both
//! faces therefore covers every block.
void FETiedFluidFSI::BuildMatrixProfile(FEGlobalMatrix& K)
{
    FEMesh& mesh = GetMesh();

    // get the DOFS
    const int dof_X = GetDOFIndex("x");
    const int dof_Y = GetDOFIndex("y");
    const int dof_Z = GetDOFIndex("z");
    const int dof_WX = GetDOFIndex("wx");
    const int dof_WY = GetDOFIndex("wy");
    const int dof_WZ = GetDOFIndex("wz");
    const int dof_EF = GetDOFIndex("ef");

    const int ndpn = 7;
    vector<int> lm(ndpn*FEElement::MAX_NODES*2);

    int npass = (m_btwo_pass ? 2 : 1);
    for (int np=0; np<npass; ++np)
    {
        FETiedFluidFSISurface& ss = (np == 0? m_ss : m_ms);

        for (int j=0; j<ss.Elements(); ++j)
        {
            FESurfaceElement& se = ss.Element(j);
            int nint = se.GaussPoints();
            int* sn = &se.m_node[0];
            for (int k=0; k<nint; ++k)
            {
                FETiedFluidFSISurface::Data& pt = static_cast<FETiedFluidFSISurface::Data&>(*se.GetMaterialPoint(k));
                FESurfaceElement* pe = pt.m_pme;
                if (pe != 0)
                {
                    FESurfaceElement& me = *pe;
                    int* mn = &me.m_node[0];

                    assign(lm, -1);

                    int nseln = se.Nodes();
                    int nmeln = me.Nodes();

                    for (int l=0; l<nseln; ++l)
                    {
                        vector<int>& id = mesh.Node(sn[l]).m_ID;
                        lm[ndpn*l  ] = id[dof_X];
                        lm[ndpn*l+1] = id[dof_Y];
                        lm[ndpn*l+2] = id[dof_Z];
                        lm[ndpn*l+3] = id[dof_WX];
                        lm[ndpn*l+4] = id[dof_WY];
                        lm[ndpn*l+5] = id[dof_WZ];
                        lm[ndpn*l+6] = id[dof_EF];
                    }

                    for (int l=0; l<nmeln; ++l)
                    {
                        vector<int>& id = mesh.Node(mn[l]).m_ID;
                        lm[ndpn*(l+nseln)  ] = id[dof_X];
                        lm[ndpn*(l+nseln)+1] = id[dof_Y];
                        lm[ndpn*(l+nseln)+2] = id[dof_Z];
                        lm[ndpn*(l+nseln)+3] = id[dof_WX];
                        lm[ndpn*(l+nseln)+4] = id[dof_WY];
                        lm[ndpn*(l+nseln)+5] = id[dof_WZ];
                        lm[ndpn*(l+nseln)+6] = id[dof_EF];
                    }

                    K.build_add(lm);
                }
            }
        }
    }

    // In an FSI-solid pair the owned fluid traction couples each FSI face's u dofs to
    // that same face's u, w and ef dofs. The tied-pair blocks above cover this only for
    // faces that actually projected onto the secondary surface, whereas the traction
    // applies on every face of the FSI surface, so add the self-coupling explicitly.
    if (m_mode == TIE_FSI_SOLID)
    {
        for (int j=0; j<m_ss.Elements(); ++j)
        {
            FESurfaceElement& se = m_ss.Element(j);
            int nseln = se.Nodes();
            int* sn = &se.m_node[0];

            assign(lm, -1);
            for (int l=0; l<nseln; ++l)
            {
                vector<int>& id = mesh.Node(sn[l]).m_ID;
                lm[ndpn*l  ] = id[dof_X];
                lm[ndpn*l+1] = id[dof_Y];
                lm[ndpn*l+2] = id[dof_Z];
                lm[ndpn*l+3] = id[dof_WX];
                lm[ndpn*l+4] = id[dof_WY];
                lm[ndpn*l+5] = id[dof_WZ];
                lm[ndpn*l+6] = id[dof_EF];
            }
            lm.resize(ndpn*nseln);
            K.build_add(lm);
            lm.resize(ndpn*FEElement::MAX_NODES*2);
        }
    }
}

//-----------------------------------------------------------------------------
//! Recalculate all three automatic penalty factors.
//!
//! NOTE: in a two-pass analysis the second pass integrates over m_ms, so the
//! integration point penalties on m_ms must be evaluated as well (otherwise they
//! retain their default value of 1).
void FETiedFluidFSI::UpdateAutoPenalty()
{
    if (m_bautopen)
    {
        CalcAutoPenalty(m_ss);
        CalcAutoViscousTractionPenalty(m_ss);
        CalcAutoNormalVelocityPenalty(m_ss);
        if (m_btwo_pass) {
            CalcAutoPenalty(m_ms);
            CalcAutoViscousTractionPenalty(m_ms);
            CalcAutoNormalVelocityPenalty(m_ms);
        }
    }
}

//-----------------------------------------------------------------------------
void FETiedFluidFSI::Activate()
{
    // don't forget to call the base class
    FEContactInterface::Activate();

    UpdateAutoPenalty();

    // project the surfaces onto each other
    // this will evaluate the gap functions in the reference configuration
    InitialProjection(m_ss, m_ms, true);
    if (m_btwo_pass) InitialProjection(m_ms, m_ss, false);

    // The owned traction load builds its element list and face orientations here.
    // Its own Activate assumes a face with one neighbour is a free surface and flips
    // the sign; a tied FSI face is exterior but is an interface with a solid, so the
    // interface orientation has to be imposed afterwards.
    if (m_pFSItrac) {
        m_pFSItrac->Activate();
        m_pFSItrac->SetInterfaceOrientation();
    }
}

//-----------------------------------------------------------------------------
//! Recompute the automatic penalties at the start of each time step, if asked.
void FETiedFluidFSI::PrepStep()
{
    if (m_bupdtpen) UpdateAutoPenalty();
}

//-----------------------------------------------------------------------------
//! The solid traction penalty eps_s has units of stress per length, and is obtained
//! from the standard elastic auto-penalty of the solid backing the surface.
void FETiedFluidFSI::CalcAutoPenalty(FETiedFluidFSISurface& s)
{
    // loop over all surface elements
    for (int i=0; i<s.Elements(); ++i)
    {
        // get the surface element
        FESurfaceElement& el = s.Element(i);

        // calculate a penalty
        double eps = AutoPenalty(el, s);

        // assign to integration points of surface element
        int nint = el.GaussPoints();
        for (int j=0; j<nint; ++j)
        {
            FETiedFluidFSISurface::Data& pt = static_cast<FETiedFluidFSISurface::Data&>(*el.GetMaterialPoint(j));
			pt.m_epss = eps;
        }
    }
}

//-----------------------------------------------------------------------------
//! Return the fluid of the fluid-FSI material backing this surface element.
//! In the FSI framework the material is a FEFluidFSI, which is not itself a
//! FEFluidMaterial; the fluid has to be extracted from it.
FEFluid* FETiedFluidFSI::GetFluid(FESurfaceElement& el)
{
    if (el.m_elem[0].pe == nullptr) return nullptr;

    // the surface may be backed by either of the two elements it separates;
    // try both before giving up
    for (int i=0; i<2; ++i)
    {
        FEElement* pe = el.m_elem[i].pe;
        if (pe == nullptr) continue;
        FEMaterial* pmat = GetFEModel()->GetMaterial(pe->GetMatID());
        FEFluidFSI* pfsi = dynamic_cast<FEFluidFSI*>(pmat);
        if (pfsi) return pfsi->Fluid();
    }

    return nullptr;
}

//-----------------------------------------------------------------------------
void FETiedFluidFSI::CalcAutoViscousTractionPenalty(FETiedFluidFSISurface& s)
{
    // loop over all surface elements
    for (int i=0; i<s.Elements(); ++i)
    {
        // get the surface element
        FESurfaceElement& el = s.Element(i);

        // calculate a penalty
        double eps = AutoViscousTractionPenalty(el, s);

        // assign to integration points of surface element
        int nint = el.GaussPoints();
        for (int j=0; j<nint; ++j)
        {
			FETiedFluidFSISurface::Data& pt = static_cast<FETiedFluidFSISurface::Data&>(*el.GetMaterialPoint(j));
			pt.m_epst = eps;
        }
    }
}

//-----------------------------------------------------------------------------
//! The viscous traction penalty has units of viscosity per length, therefore it is
//! scaled by the ratio of the fluid viscosity to the element thickness.
double FETiedFluidFSI::AutoViscousTractionPenalty(FESurfaceElement& el, FETiedFluidFSISurface& s)
{
    // get the solid element attached to the surface element
    if (el.m_elem[0].pe == nullptr) return 1.0;
    FESolidElement& sel = static_cast<FESolidElement&>(*el.m_elem[0].pe);

    // get the fluid of the FSI material for that element
    FEFluid* pfluid = GetFluid(el);
    if (pfluid == nullptr) return 1.0;

    // get the viscous part of this fluid
    FEViscousFluid* pvfluid = pfluid->GetViscous();
    if (pvfluid == nullptr) return 1.0;

    // evaluate the viscosity for each material point and get its average
    double eta = 0;
    int nint = sel.GaussPoints();
    for (int i=0; i<nint; ++i) {
        FEMaterialPoint* mp = sel.GetMaterialPoint(i);
        eta += pvfluid->ShearViscosity(*mp);
    }
    eta /= nint;

    // get the element thickness (both measures evaluated in the reference configuration)
    double area = s.GetArea(el, true);
    double vol = s.GetVolume(sel);
    if (area <= 0) return 1.0;
    double h = vol/area;
    if (h <= 0) return 1.0;

    return eta/h;
}

//-----------------------------------------------------------------------------
void FETiedFluidFSI::CalcAutoNormalVelocityPenalty(FETiedFluidFSISurface& s)
{
    // loop over all surface elements
    for (int i=0; i<s.Elements(); ++i)
    {
        // get the surface element
        FESurfaceElement& el = s.Element(i);

        // calculate a penalty
        double eps = AutoNormalVelocityPenalty(el, s);

        // assign to integration points of surface element
        int nint = el.GaussPoints();
        for (int j=0; j<nint; ++j)
        {
            FETiedFluidFSISurface::Data& pt = static_cast<FETiedFluidFSISurface::Data&>(*el.GetMaterialPoint(j));
            pt.m_epsn = eps;
        }
    }
}

//-----------------------------------------------------------------------------
//! The normal velocity penalty has units of velocity. Without augmentation the
//! converged dilatation gap is pi = w_n/eps_n, so eps_n must be large compared to the
//! throughflow velocity divided by the acceptable gap. Choosing eps_n = h/tau = h*K/eta
//! makes the pressure drop across the interface eta*w/h, i.e. a small fraction of the
//! channel pressure drop, so the tie is transparent. For an inviscid fluid there is no
//! viscous time scale and we fall back on the acoustic impedance sqrt(K/rho).
double FETiedFluidFSI::AutoNormalVelocityPenalty(FESurfaceElement& el, FETiedFluidFSISurface& s)
{
    // get the solid element attached to the surface element
    if (el.m_elem[0].pe == nullptr) return 1.0;
    FESolidElement& sel = static_cast<FESolidElement&>(*el.m_elem[0].pe);

    // get the fluid of the FSI material for that element
    FEFluid* pfluid = GetFluid(el);
    if (pfluid == nullptr) return 1.0;

    // bulk modulus
    double K = pfluid->m_k;
    if (K <= 0) return 1.0;

    // evaluate the viscosity for each material point and get its average
    FEViscousFluid* pvfluid = pfluid->GetViscous();
    double eta = 0;
    if (pvfluid) {
        int nint = sel.GaussPoints();
        for (int i=0; i<nint; ++i) {
            FEMaterialPoint* mp = sel.GetMaterialPoint(i);
            eta += pvfluid->ShearViscosity(*mp);
        }
        eta /= nint;
    }

    // get the element thickness (both measures evaluated in the reference configuration)
    double area = s.GetArea(el, true);
    double vol = s.GetVolume(sel);
    double h = (area > 0 ? vol/area : 0);

    if ((eta > 0) && (h > 0)) {
        // viscous relaxation scaling h/tau, with tau = eta/K
        double tau = eta/K;
        return h/tau;
    }

    // inviscid fluid: fall back on the acoustic impedance of the fluid
    double rho = pfluid->m_rhor;
    if (rho > 0) return sqrt(K/rho);

    return 1.0;
}

//-----------------------------------------------------------------------------
// Perform initial projection between tied surfaces in reference configuration
void FETiedFluidFSI::InitialProjection(FETiedFluidFSISurface& ss, FETiedFluidFSISurface& ms, bool bfirst)
{
    FESurfaceElement* pme;
    vec3d r, nu;
    double rs[2];

    // initialize projection data
    FENormalProjection np(ms);
    np.SetTolerance(m_stol);
    np.SetSearchRadius(m_srad);
    np.Init();

    // projection diagnostics
    int nproj = 0, nfail = 0;
    double maxgap = 0;

    // loop over all integration points
    for (int i=0; i<ss.Elements(); ++i)
    {
        FESurfaceElement& el = ss.Element(i);

        int nint = el.GaussPoints();

        for (int j=0; j<nint; ++j)
        {
            // calculate the global position of the integration point
            r = ss.Local2Global(el, j);

            // calculate the normal at this integration point
            nu = ss.SurfaceNormal(el, j);

            // find the intersection point with the secondary surface
            pme = np.Project2(r, nu, rs);

            FETiedFluidFSISurface::Data& pt = static_cast<FETiedFluidFSISurface::Data&>(*el.GetMaterialPoint(j));
			pt.m_pme = pme;
            pt.m_nu = nu;
            pt.m_rs[0] = rs[0];
            pt.m_rs[1] = rs[1];
            if (pme)
            {
                // the node could potentially be tied
                // find the global location of the intersection point
                vec3d q = ms.Local2Global(*pme, rs[0], rs[1]);

                // calculate the gap function
                pt.m_Gap = q - r;
                ++nproj;
                maxgap = max(maxgap, pt.m_Gap.norm());
            }
            else
            {
                // the integration point could not be projected onto the opposing surface
                pt.m_Gap = vec3d(0,0,0);
                ++nfail;
            }
        }
    }

    // report the outcome of the projection. Integration points that could not be
    // projected are left untied: they transmit no solid traction, no viscous traction
    // and no normal flow, i.e. they revert to a free, impermeable, frictionless wall.
    feLog(" tied fluid-FSI interface # %d: %s pass\n", GetID(), (bfirst ? "primary" : "secondary"));
    feLog("    tied integration points  : %d\n", nproj);
    feLog("    maximum initial gap      : %15le\n", maxgap);
    if (nproj == 0)
    {
        std::string name = GetName();
        feLogWarning("No contact pairs found for tied interface \"%s\".\nThis contact interface may not have any effect.", name.c_str());
    }
    else if (nfail > 0) {
        feLogWarning("Tied fluid-FSI interface %d: %d integration point(s) could not be projected onto the\n"
                     "opposing surface. These points remain untied and transmit no traction or flow.\n"
                     "Consider increasing search_radius or search_tol.", GetID(), nfail);
    }
}

//-----------------------------------------------------------------------------
//! Evaluate the three gap functions and the corresponding tractions.
//!
//! All gaps are oriented as (secondary - primary):
//!   g   = (q - r) - Gap    solid displacement gap, eq. (7.9.5)
//!   g^w = w(2) - w(1)      relative fluid velocity gap, eq. (7.9.4)
//!   pi  = J(2) - J(1)      dilatation gap
//!
//! and the tractions follow the augmented Lagrangian forms
//!   t^s = lambda_s + eps_s*g, t^tau = lambda_t + eps_t*g^w, w_n = lambda_p + eps_n*pi,
//! which reduce to the penalty forms when the multipliers are zero. Against a solid
//! wall t^tau is additionally projected onto the normal, see below.
//!
//! In TIE_FSI_SOLID mode the secondary surface is a plain solid, which has neither w
//! nor J. It is treated as a wall with w(2) = 0 and only the normal part of the
//! resulting gap is kept, g^w = -(w(1).n) n, which enforces impermeability w.n = 0
//! and leaves tangential slip free; and the dilatation gap is dropped altogether,
//! since J is unconstrained at an FSI-solid boundary.
void FETiedFluidFSI::ProjectSurface(FETiedFluidFSISurface& ss, FETiedFluidFSISurface& ms)
{
    FEMesh& mesh = GetMesh();
    FESurfaceElement* pme;
    vec3d r;
    double alpha = GetFEModel()->GetTime().alphaf;

    vec3d  wt[FEElement::MAX_NODES], wp[FEElement::MAX_NODES], w1;
    double et[FEElement::MAX_NODES], ep[FEElement::MAX_NODES], e1;

    // loop over all integration points
    for (int i=0; i<ss.Elements(); ++i)
    {
        FESurfaceElement& el = ss.Element(i);

        int ne = el.Nodes();
        int nint = el.GaussPoints();

        // get the nodal relative fluid velocities and dilatations
        for (int j=0; j<ne; ++j) {
            FENode& node = mesh.Node(el.m_node[j]);
            wt[j] = node.get_vec3d(m_dofW[0], m_dofW[1], m_dofW[2]);
            wp[j] = node.get_vec3d_prev(m_dofW[0], m_dofW[1], m_dofW[2]);
            et[j] = node.get(m_dofEF);
            ep[j] = node.get_prev(m_dofEF);
        }

        for (int j=0; j<nint; ++j)
        {
            FETiedFluidFSISurface::Data& pt = static_cast<FETiedFluidFSISurface::Data&>(*el.GetMaterialPoint(j));

            // Calculate the global position of the integration point.
            // NOTE: the position, the normal, and (in LoadVector and StiffnessMatrix)
            // the covariant basis vectors and J_eta, are all evaluated in the
            // alpha_f-weighted intermediate configuration, because that is the
            // configuration about which the alpha_f prefactors of eqs. (7.9.13)-(7.9.16)
            // linearize. Evaluating them at the end of the step instead would leave the
            // u-columns of the tangent short by a factor alpha_f whenever alpha_f < 1
            // (i.e. whenever numerical damping is requested through rhoi > 0).
            r = ss.Local2GlobalAlpha(el, j, alpha);

            // calculate the normal at this integration point
            pt.m_nu = ss.SurfaceNormalAlpha(el, j, alpha);

            // relative fluid velocity and dilatation at the integration point
            w1 = el.eval(wt, j)*alpha + el.eval(wp, j)*(1-alpha);
            e1 = el.eval(et, j)*alpha + el.eval(ep, j)*(1-alpha);

            // if this node is tied, evaluate gap functions
            pme = pt.m_pme;
            if (pme)
            {
                // find the global location of the intersection point
                vec3d q = ms.Local2GlobalAlpha(*pme, pt.m_rs[0], pt.m_rs[1], alpha);

                // solid displacement gap, measured from the initial gap
                pt.m_dg = (q - r) - pt.m_Gap;
				pt.m_gap = pt.m_dg.norm();

                // relative fluid velocity and dilatation gaps
                if (ms.IsFSI())
                {
                    vec3d wmt[FEElement::MAX_NODES], wmp[FEElement::MAX_NODES];
                    double emt[FEElement::MAX_NODES], emp[FEElement::MAX_NODES];
                    for (int k=0; k<pme->Nodes(); ++k) {
                        FENode& node = mesh.Node(pme->m_node[k]);
                        wmt[k] = node.get_vec3d(m_dofW[0], m_dofW[1], m_dofW[2]);
                        wmp[k] = node.get_vec3d_prev(m_dofW[0], m_dofW[1], m_dofW[2]);
                        emt[k] = node.get(m_dofEF);
                        emp[k] = node.get_prev(m_dofEF);
                    }
                    vec3d w2 = pme->eval(wmt, pt.m_rs[0], pt.m_rs[1])*alpha + pme->eval(wmp, pt.m_rs[0], pt.m_rs[1])*(1-alpha);
                    pt.m_gw = w2 - w1;

                    double e2 = pme->eval(emt, pt.m_rs[0], pt.m_rs[1])*alpha + pme->eval(emp, pt.m_rs[0], pt.m_rs[1])*(1-alpha);
                    pt.m_Jg = e2 - e1;
                }
                else
                {
                    // Solid wall: w(2) = 0, so the gap would be -w(1); keeping only its
                    // normal part enforces impermeability w.n = 0 and leaves tangential
                    // slip free. J is unconstrained here.
                    pt.m_gw = pt.m_nu*(-(w1*pt.m_nu));
                    pt.m_Jg = 0.0;
                }

                // penalty factors
                double epss = m_epss*pt.m_epss;
                double epst = m_epst*pt.m_epst;
                double epsn = m_epsn*pt.m_epsn;

                // remember w on this surface; the tangent of the wall constraint needs
                // the full vector, not just the normal component retained in m_gw
                pt.m_wr = w1;

                // tractions and normal velocity (augmented Lagrangian form)
                pt.m_ts = pt.m_Lmd + pt.m_dg*epss;
                if (ms.IsFSI())
                    pt.m_tv = pt.m_Lmt + pt.m_gw*epst;
                else
                    // purely normal: project the multiplier too, so that a normal that
                    // has rotated since the last augmentation cannot leave a tangential
                    // residue in the traction
                    pt.m_tv = pt.m_nu*((pt.m_Lmt*pt.m_nu) + (pt.m_gw*pt.m_nu)*epst);
                pt.m_wn = (ms.IsFSI() ? pt.m_Lmp + pt.m_Jg*epsn : 0.0);
            }
            else
            {
                // the node is not tied
                pt.m_dg = vec3d(0,0,0);
				pt.m_gap = 0.0;
                pt.m_gw = vec3d(0,0,0);
                pt.m_wr = vec3d(0,0,0);
                pt.m_Jg = 0.0;
                pt.m_ts = vec3d(0,0,0);
                pt.m_tv = vec3d(0,0,0);
                pt.m_wn = 0.0;
            }
        }
    }
}

//-----------------------------------------------------------------------------

void FETiedFluidFSI::Update()
{
    // project the surfaces onto each other
    // this will update the gap functions as well
    ProjectSurface(m_ss, m_ms);
    if (m_btwo_pass) ProjectSurface(m_ms, m_ss);
}

//-----------------------------------------------------------------------------
//! Contact forces, eq. (7.9.10)-(7.9.11).
//!
//! Seven dofs per node, ordered u (3), w (3), J (1):
//!     f^s_a = N^(1)_a t^s      f^s_b = -N^(2)_b t^s
//!     f_a   = N^(1)_a t^tau    f_b   = -N^(2)_b t^tau
//!     f^J_a = N^(1)_a w_n      f^J_b = -N^(2)_b w_n
//! each weighted by W_k J_eta.
//!
//! Against a solid secondary surface only the three u rows exist on that side, and
//! w_n is zero, so the secondary contributes f^s_b alone and the primary J row
//! vanishes. The fluid traction is added at the end by the owned traction load.
void FETiedFluidFSI::LoadVector(FEGlobalVector& R, const FETimeInfo& tp)
{
    vector<int> LM1, LM2, LM, en;
    vector<double> fe;
    const int MI = FEElement::MAX_INTPOINTS;
    const int MN = FEElement::MAX_NODES;
    double detJ[MI], w[MI], *H1, H2[MN];
    vec3d s1v[MN], s2v[MN], f1[MN], f2[MN];
    double w1[MN], w2[MN];

    // zero the nodal forces on both surfaces
    zero(m_ss.m_Fn);
    zero(m_ms.m_Fn);

    // loop over the nr of passes
    int npass = (m_btwo_pass?2:1);
    for (int np=0; np<npass; ++np)
    {
        // get primary and secondary surfaces
        FETiedFluidFSISurface& ss = (np == 0? m_ss : m_ms);
        FETiedFluidFSISurface& ms = (np == 0? m_ms : m_ss);

        // dofs per node on each side: 7 for an FSI surface, 3 for a plain solid
        const int ndpn1 = ss.DofsPerNode();
        const int ndpn2 = ms.DofsPerNode();

        // loop over all elements of primary surface
        for (int i=0; i<ss.Elements(); ++i)
        {
            // get the surface element
            FESurfaceElement& se1 = ss.Element(i);

            // get the nr of nodes and integration points
            int neln1 = se1.Nodes();
            int nint1 = se1.GaussPoints();

            // copy the LM vector; we'll need it later
            ss.UnpackLM(se1, LM1);

            // we calculate all the metrics we need before we
            // calculate the nodal forces
            for (int j=0; j<nint1; ++j)
            {
                // get the base vectors in the intermediate configuration
                vec3d g[2];
                ss.CoBaseVectorsAlpha(se1, j, tp.alphaf, g);

                // jacobians: J_eta = |g1xg2|
                detJ[j] = (g[0] ^ g[1]).norm();

                // integration weights
                w[j] = se1.GaussWeights()[j];
            }

            // loop over all integration points
            // note that we are integrating over the current surface
            for (int j=0; j<nint1; ++j)
            {
				FETiedFluidFSISurface::Data& pt = static_cast<FETiedFluidFSISurface::Data&>(*se1.GetMaterialPoint(j));

                // get the secondary surface element
                FESurfaceElement* pme = pt.m_pme;
                if (pme)
                {
                    // get the secondary surface element
                    FESurfaceElement& se2 = *pme;

                    // get the nr of secondary element nodes
                    int neln2 = se2.Nodes();

                    // copy LM vector
                    ms.UnpackLM(se2, LM2);

                    // calculate degrees of freedom
                    int ndof = ndpn1*neln1 + ndpn2*neln2;
                    int off2 = ndpn1*neln1;

                    // build the LM vector
                    LM.resize(ndof);
                    for (int a=0; a<neln1; ++a)
                        for (int k=0; k<ndpn1; ++k) LM[ndpn1*a+k] = LM1[ndpn1*a+k];

                    for (int b=0; b<neln2; ++b)
                        for (int k=0; k<ndpn2; ++k) LM[off2+ndpn2*b+k] = LM2[ndpn2*b+k];

                    // build the en vector
                    en.resize(neln1+neln2);
                    for (int a=0; a<neln1; ++a) en[a      ] = se1.m_node[a];
                    for (int b=0; b<neln2; ++b) en[b+neln1] = se2.m_node[b];

                    // get primary element shape functions
                    H1 = se1.H(j);

                    // get secondary element shape functions
                    double r = pt.m_rs[0];
                    double s = pt.m_rs[1];
                    se2.shape_fnc(H2, r, s);

                    // tractions and normal velocity
                    // (evaluated in ProjectSurface, called from Update)
                    vec3d ts = pt.m_ts;
                    vec3d tv = pt.m_tv;
                    double wn = pt.m_wn;

                    // calculate the force vector
                    fe.resize(ndof);
                    zero(fe);

                    for (int a=0; a<neln1; ++a) {
                        s1v[a] = ts*H1[a];
                        f1[a]  = tv*H1[a];
                        w1[a]  = wn*H1[a];
                    }
                    for (int b=0; b<neln2; ++b) {
                        s2v[b] = -ts*H2[b];
                        f2[b]  = -tv*H2[b];
                        w2[b]  = -wn*H2[b];
                    }

                    for (int a=0; a<neln1; ++a)
                    {
                        fe[ndpn1*a    ] += s1v[a].x*detJ[j]*w[j];
                        fe[ndpn1*a + 1] += s1v[a].y*detJ[j]*w[j];
                        fe[ndpn1*a + 2] += s1v[a].z*detJ[j]*w[j];
                        fe[ndpn1*a + 3] += f1[a].x*detJ[j]*w[j];
                        fe[ndpn1*a + 4] += f1[a].y*detJ[j]*w[j];
                        fe[ndpn1*a + 5] += f1[a].z*detJ[j]*w[j];
                        // against a solid wall w_n is zero, so this row vanishes anyway
                        fe[ndpn1*a + 6] += w1[a]*detJ[j]*w[j];

                        // accumulate the nodal solid force for reporting
                        ss.m_Fn[se1.m_lnode[a]] += s1v[a]*(detJ[j]*w[j]);
                    }
                    for (int b=0; b<neln2; ++b) {
                        fe[off2 + ndpn2*b    ] += s2v[b].x*detJ[j]*w[j];
                        fe[off2 + ndpn2*b + 1] += s2v[b].y*detJ[j]*w[j];
                        fe[off2 + ndpn2*b + 2] += s2v[b].z*detJ[j]*w[j];
                        if (ms.IsFSI()) {
                            // a solid secondary surface has no w and no J rows
                            fe[off2 + ndpn2*b + 3] += f2[b].x*detJ[j]*w[j];
                            fe[off2 + ndpn2*b + 4] += f2[b].y*detJ[j]*w[j];
                            fe[off2 + ndpn2*b + 5] += f2[b].z*detJ[j]*w[j];
                            fe[off2 + ndpn2*b + 6] += w2[b]*detJ[j]*w[j];
                        }

                        ms.m_Fn[se2.m_lnode[b]] += s2v[b]*(detJ[j]*w[j]);
                    }

                    // assemble the global residual
                    R.Assemble(en, LM, fe);
                }
            }
        }
    }

    // Transfer the fluid traction to the solid. The tie above only balances t^s(1)
    // against t^solid(2); sigma^f.n enters neither u-equation, so without this the
    // fluid pressure and viscous stress would have no effect on the interface. The
    // load acts on the FSI surface's own u dofs and the displacement tie carries it
    // across, reproducing the conforming-mesh balance F^s(1) + F^solid(2) = T^f.
    if (m_pFSItrac) m_pFSItrac->LoadVector(R);
}

//-----------------------------------------------------------------------------
//! Contact stiffness, eqs. (7.9.13)-(7.9.17).
//!
//! The block structure of eq. (7.9.12), with rows [dv^s_a, dw_a, dJ_a] and columns
//! [Du_c, Dw_c, DJ_c], is
//!
//!     (1,1) and (2,1):  [ K^s    0    0 ]        (1,2) and (2,2):  [ K^s  0    0 ]
//!                       [ K^fs   K    0 ]                          [ 0    K    0 ]
//!                       [ k^fs   0    k ]                          [ 0    0    k ]
//!
//! The K^fs and k^fs blocks appear only in the columns of Du^(1), because only the
//! primary surface geometry (through J_eta and n^(1)) varies with the solid motion.
//! Those blocks are also what makes this tangent non-symmetric.
void FETiedFluidFSI::StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp)
{
    vector<int> LM1, LM2, LM, en;
    const int MI = FEElement::MAX_INTPOINTS;
    const int MN = FEElement::MAX_NODES;
    double detJ[MI], w[MI], *H1, H2[MN];
    vec3d gtan[MI][2];
    FEElementMatrix ke;

    // scratch for the dn/du term of the wall constraint, filled per integration point
    // in TIE_FSI_SOLID mode only. Declared here rather than in the integration point
    // loop so that mat3d's zero-initialising default constructor does not run over the
    // whole array at every point, including in FSI-FSI mode where it is never used.
    mat3d dtau[MN];

    double alpha = tp.alphaf;

    // Set the higher order stiffness multiplier. Following the convention of the other
    // FEBio tied interfaces, a negative knmult = -n means the higher order (geometric)
    // terms are switched on only after n stiffness reformations, which is useful for
    // comparing the contribution of the different stiffness terms.
    int nref = GetSolver()->m_nref;
    double knmult = m_knmult;
    if (m_knmult < 0)
    {
        int ni = int(-m_knmult);
        if (nref >= ni)
        {
            knmult = 1;
            feLog("Higher order stiffness terms included.\n");
        }
        else knmult = 0;
    }

    // do single- or two-pass
    int npass = (m_btwo_pass?2:1);
    for (int np=0; np < npass; ++np)
    {
		// get primary and secondary surfaces
		FETiedFluidFSISurface& ss = (np == 0? m_ss : m_ms);
        FETiedFluidFSISurface& ms = (np == 0? m_ms : m_ss);

        // dofs per node on each side: 7 for an FSI surface, 3 for a plain solid.
        // When the secondary side is a solid it has no w and no J rows, so the blocks
        // K^(2,1), k^(2,1), K^(1,2), K^(2,2), k^(1,2), k^(2,2) and the whole of
        // K^fs(2,1)/k^fs(2,1) have no rows or columns to occupy and are skipped.
        const int ndpn1 = ss.DofsPerNode();
        const int ndpn2 = ms.DofsPerNode();
        const bool bmsFSI = ms.IsFSI();

        // loop over all elements of primary surface
        for (int i=0; i<ss.Elements(); ++i)
        {
            // get the next element
            FESurfaceElement& se1 = ss.Element(i);

            // get nr of nodes and integration points
            int neln1 = se1.Nodes();
            int nint1 = se1.GaussPoints();

            // copy the LM vector
            ss.UnpackLM(se1, LM1);

            // we calculate all the metrics we need before we
            // calculate the nodal forces
            for (int j=0; j<nint1; ++j)
            {
                // get the base vectors in the intermediate configuration, the same
                // configuration ProjectSurface and LoadVector use
                ss.CoBaseVectorsAlpha(se1, j, alpha, gtan[j]);

                // jacobians: J_eta = |g1xg2|
                detJ[j] = (gtan[j][0] ^ gtan[j][1]).norm();

                // integration weights
                w[j] = se1.GaussWeights()[j];
            }

            // loop over all integration points
            for (int j=0; j<nint1; ++j)
            {
                FETiedFluidFSISurface::Data& pt = static_cast<FETiedFluidFSISurface::Data&>(*se1.GetMaterialPoint(j));

                // get the secondary element
                FESurfaceElement* pme = pt.m_pme;
                if (pme)
                {
                    FESurfaceElement& se2 = *pme;

                    // get the nr of secondary nodes
                    int neln2 = se2.Nodes();

                    // copy the LM vector
                    ms.UnpackLM(se2, LM2);

                    int ndof = ndpn1*neln1 + ndpn2*neln2;
                    int off2 = ndpn1*neln1;

                    // build the LM vector
                    LM.resize(ndof);

                    for (int a=0; a<neln1; ++a)
                        for (int k=0; k<ndpn1; ++k) LM[ndpn1*a+k] = LM1[ndpn1*a+k];

                    for (int b=0; b<neln2; ++b)
                        for (int k=0; k<ndpn2; ++k) LM[off2+ndpn2*b+k] = LM2[ndpn2*b+k];

                    // build the en vector
                    en.resize(neln1+neln2);
                    for (int a=0; a<neln1; ++a) en[a      ] = se1.m_node[a];
                    for (int b=0; b<neln2; ++b) en[b+neln1] = se2.m_node[b];

                    // primary shape functions and their derivatives
                    H1 = se1.H(j);
                    double* Gr1 = se1.Gr(j);
                    double* Gs1 = se1.Gs(j);

                    // secondary shape functions
                    double r = pt.m_rs[0];
                    double s = pt.m_rs[1];
                    se2.shape_fnc(H2, r, s);

                    // penalties
                    double epss = m_epss*pt.m_epss;
                    double epst = m_epst*pt.m_epst;
                    double epsn = m_epsn*pt.m_epsn;

                    // tractions, normal velocity and normal
                    vec3d ts = pt.m_ts;
                    vec3d tv = pt.m_tv;
                    double wn = pt.m_wn;
                    vec3d nu = pt.m_nu;

                    // the geometric terms of eq. (7.9.14) all involve
                    //     A^(1)_c . n^(1),  with  A^(1)_c = A{ dN_c/deta2 g1 - dN_c/deta1 g2 }
                    // the antisymmetric tensor whose dual vector is in the braces, eq. (7.9.17).
                    // Since A(a).b = a x b, we form the cross product directly.
                    vec3d Acn[MN];
                    vec3d Acv[MN];
                    for (int c=0; c<neln1; ++c) {
                        Acv[c] = gtan[j][0]*Gs1[c] - gtan[j][1]*Gr1[c];
                        Acn[c] = Acv[c] ^ nu;
                    }

                    // The A_c terms are the non-symmetric part of the tangent, scaled by
                    // knmult. Under symmetric_stiffness the solid block keeps the
                    // symmetric part of its A_c term, as the other FEBio tied interfaces
                    // do, whereas the K^fs and k^fs blocks are dropped: those couple
                    // different variable groups (w-rows and J-rows against u-columns) and
                    // have no diagonal counterpart to symmetrize against.
                    double sscale = knmult/detJ[j];
                    double gscale = (m_bsymm ? 0.0 : sscale);

                    // Against a solid wall t^tau = (lambda.n + eps_t*g_n) n depends on u
                    // through the normal as well as through J_eta. With
                    //     dn/du_c = alpha (A_c + n (x) (A_c.n)) / J_eta = alpha P_c ,
                    //     g_n = -w.n ,   t_n = t^tau . n ,
                    // the extra contribution to the w-row is
                    //     dt^tau/du_c = n (x) (P_c^T (lambda - eps_t w)) + t_n P_c ,
                    // which depends only on c, so evaluate it once per c here.
                    if (bmsFSI == false)
                    {
                        double tn = tv*nu;
                        for (int c=0; c<neln1; ++c) {
                            mat3d Pc = (mat3d(mat3da(Acv[c])) + (nu & Acn[c]))*gscale;
                            dtau[c] = (nu & (Pc.transpose()*(pt.m_Lmt - pt.m_wr*epst))) + Pc*tn;
                        }
                    }

                    // create the stiffness matrix
                    ke.resize(ndof, ndof); ke.zero();

                    double Jw = detJ[j]*w[j];

                    //------------------------------------
                    // rows on the primary surface
                    for (int a=0; a<neln1; ++a) {
                        for (int c=0; c<neln1; ++c)
                        {
                            // eq. (7.9.13)_1 and (7.9.14)_1,3
                            mat3d Ats11 = (ts & Acn[c])*sscale;
                            if (m_bsymm) Ats11 = (Ats11 + Ats11.transpose())*0.5;
                            mat3d Ks11 = -(Ats11 + mat3dd(epss*H1[c]))*(alpha*H1[a]*Jw);
                            mat3d Kfs11 = (tv & Acn[c])*(-alpha*H1[a]*gscale*Jw);
                            vec3d kfs11 = Acn[c]*(-alpha*wn*gscale*H1[a]*Jw);
                            double k11 = -alpha*epsn*H1[a]*H1[c]*Jw;

                            // The w-w block. Against an FSI surface the constraint is
                            // [[w]] = 0 and the block is isotropic, eq. (7.9.15)_1.
                            // Against a solid wall only w.n is constrained, so the
                            // traction is t^tau = eps_t (-w.n) n and its derivative
                            // carries the projector n (x) n in place of the identity.
                            mat3d K11 = (bmsFSI ? mat3d(mat3dd(1.0)) : (nu & nu));
                            K11 *= -alpha*epst*H1[a]*H1[c]*Jw;

                            // add the dn/du contribution of the wall constraint. NOTE the
                            // sign: dtau already is dt^tau/du_c, so the total w-row
                            // derivative is W N_a [J_eta dt^tau/du_c - t^tau (x) A_c n],
                            // and this term enters with a PLUS. (Ks11 above carries a
                            // leading minus only because dt^s/du_c = -alpha eps_s N_c I.)
                            if (bmsFSI == false) Kfs11 += dtau[c]*(alpha*H1[a]*Jw);

                            for (int k=0; k<3; ++k) {
                                for (int l=0; l<3; ++l) {
                                    ke[ndpn1*a + k    ][ndpn1*c + l] -= Ks11(k,l);
                                    ke[ndpn1*a + k + 3][ndpn1*c + l] -= Kfs11(k,l);
                                    ke[ndpn1*a + k + 3][ndpn1*c + l + 3] -= K11(k,l);
                                }
                            }

                            // the J row exists only when the dilatation is constrained,
                            // i.e. when the opposing surface carries a dilatation too
                            if (bmsFSI) {
                                ke[ndpn1*a + 6][ndpn1*c    ] -= kfs11.x;
                                ke[ndpn1*a + 6][ndpn1*c + 1] -= kfs11.y;
                                ke[ndpn1*a + 6][ndpn1*c + 2] -= kfs11.z;
                                ke[ndpn1*a + 6][ndpn1*c + 6] -= k11;
                            }
                        }
                        for (int d=0; d<neln2; ++d)
                        {
                            // eq. (7.9.13)_2, (7.9.15)_2 and (7.9.16)_2
                            mat3dd Ks12(alpha*epss*H1[a]*H2[d]*Jw);

                            ke[ndpn1*a    ][off2 + ndpn2*d    ] -= Ks12.xx();
                            ke[ndpn1*a + 1][off2 + ndpn2*d + 1] -= Ks12.yy();
                            ke[ndpn1*a + 2][off2 + ndpn2*d + 2] -= Ks12.zz();

                            // a solid secondary surface has no w and no J columns
                            if (bmsFSI) {
                                mat3dd K12(alpha*epst*H1[a]*H2[d]*Jw);
                                double k12 = alpha*epsn*H1[a]*H2[d]*Jw;

                                ke[ndpn1*a + 3][off2 + ndpn2*d + 3] -= K12.xx();
                                ke[ndpn1*a + 4][off2 + ndpn2*d + 4] -= K12.yy();
                                ke[ndpn1*a + 5][off2 + ndpn2*d + 5] -= K12.zz();
                                ke[ndpn1*a + 6][off2 + ndpn2*d + 6] -= k12;
                            }
                        }
                    }

                    //------------------------------------
                    // rows on the secondary surface
                    for (int b=0; b<neln2; ++b) {
                        for (int c=0; c<neln1; ++c)
                        {
                            // eq. (7.9.13)_3 and (7.9.14)_2,4
                            mat3d Ats21 = (ts & Acn[c])*sscale;
                            if (m_bsymm) Ats21 = (Ats21 + Ats21.transpose())*0.5;
                            mat3d Ks21 = (Ats21 + mat3dd(epss*H1[c]))*(alpha*H2[b]*Jw);

                            for (int k=0; k<3; ++k)
                                for (int l=0; l<3; ++l)
                                    ke[off2 + ndpn2*b + k][ndpn1*c + l] -= Ks21(k,l);

                            // a solid secondary surface has no w and no J rows
                            if (bmsFSI) {
                                mat3d Kfs21 = (tv & Acn[c])*(alpha*H2[b]*gscale*Jw);
                                vec3d kfs21 = Acn[c]*(alpha*wn*gscale*H2[b]*Jw);
                                mat3dd K21(alpha*epst*H2[b]*H1[c]*Jw);
                                double k21 = alpha*epsn*H2[b]*H1[c]*Jw;

                                for (int k=0; k<3; ++k)
                                    for (int l=0; l<3; ++l)
                                        ke[off2 + ndpn2*b + k + 3][ndpn1*c + l] -= Kfs21(k,l);

                                ke[off2 + ndpn2*b + 6][ndpn1*c    ] -= kfs21.x;
                                ke[off2 + ndpn2*b + 6][ndpn1*c + 1] -= kfs21.y;
                                ke[off2 + ndpn2*b + 6][ndpn1*c + 2] -= kfs21.z;

                                ke[off2 + ndpn2*b + 3][ndpn1*c + 3] -= K21.xx();
                                ke[off2 + ndpn2*b + 4][ndpn1*c + 4] -= K21.yy();
                                ke[off2 + ndpn2*b + 5][ndpn1*c + 5] -= K21.zz();
                                ke[off2 + ndpn2*b + 6][ndpn1*c + 6] -= k21;
                            }
                        }
                        for (int d=0; d<neln2; ++d)
                        {
                            // eq. (7.9.13)_4, (7.9.15)_4 and (7.9.16)_4
                            mat3dd Ks22(-alpha*epss*H2[b]*H2[d]*Jw);

                            ke[off2 + ndpn2*b    ][off2 + ndpn2*d    ] -= Ks22.xx();
                            ke[off2 + ndpn2*b + 1][off2 + ndpn2*d + 1] -= Ks22.yy();
                            ke[off2 + ndpn2*b + 2][off2 + ndpn2*d + 2] -= Ks22.zz();

                            if (bmsFSI) {
                                mat3dd K22(-alpha*epst*H2[b]*H2[d]*Jw);
                                double k22 = -alpha*epsn*H2[b]*H2[d]*Jw;

                                ke[off2 + ndpn2*b + 3][off2 + ndpn2*d + 3] -= K22.xx();
                                ke[off2 + ndpn2*b + 4][off2 + ndpn2*d + 4] -= K22.yy();
                                ke[off2 + ndpn2*b + 5][off2 + ndpn2*d + 5] -= K22.zz();
                                ke[off2 + ndpn2*b + 6][off2 + ndpn2*d + 6] -= k22;
                            }
                        }
                    }

                    // assemble the global stiffness
                    ke.SetNodes(en);
                    ke.SetIndices(LM);
                    LS.Assemble(ke);
                }
            }
        }
    }

    // tangent of the fluid traction transfer (couples the FSI surface's u rows to its
    // own u, w and ef columns)
    if (m_pFSItrac) m_pFSItrac->StiffnessMatrix(LS);
}

//-----------------------------------------------------------------------------
//! Augmented Lagrangian update.
//!
//! The solid multiplier follows eq. (7.9.7), lambda_s <- lambda_s + eps_s*g. The
//! viscous traction and normal velocity multipliers are the natural extension of
//! the tied fluid interface of Section 7.8, and are updated the same way.
bool FETiedFluidFSI::Augment(int naug, const FETimeInfo& tp)
{
    // make sure we need to augment
	if (m_laugon != FECore::AUGLAG_METHOD) return true;

    bool bconv = true;

    int N1 = m_ss.Elements();
    int N2 = m_ms.Elements();

    // --- c a l c u l a t e   i n i t i a l   n o r m s ---
    double normS0 = 0, normL0 = 0, normJ0 = 0;
    for (int i=0; i<N1; ++i)
    {
		FESurfaceElement& s1 = m_ss.Element(i);
        for (int j=0; j<s1.GaussPoints(); ++j)
        {
			FETiedFluidFSISurface::Data& d1 = static_cast<FETiedFluidFSISurface::Data&>(*s1.GetMaterialPoint(j));
			normS0 += d1.m_Lmd*d1.m_Lmd;
			normL0 += d1.m_Lmt*d1.m_Lmt;
            normJ0 += d1.m_Lmp*d1.m_Lmp;
        }
    }
    for (int i=0; i<N2; ++i)
    {
		FESurfaceElement& s2 = m_ms.Element(i);
        for (int j=0; j<s2.GaussPoints(); ++j)
        {
			FETiedFluidFSISurface::Data& d2 = static_cast<FETiedFluidFSISurface::Data&>(*s2.GetMaterialPoint(j));
			normS0 += d2.m_Lmd*d2.m_Lmd;
			normL0 += d2.m_Lmt*d2.m_Lmt;
            normJ0 += d2.m_Lmp*d2.m_Lmp;
        }
    }
    normS0 = sqrt(normS0);
    normL0 = sqrt(normL0);
    normJ0 = sqrt(normJ0);

    // b. gap components
    // (are calculated during update)
    double maxdg = 0;
    double maxgap = 0;
    double maxJg = 0;

    // is the dilatation constrained? Only when both surfaces carry a dilatation.
    const bool bJ = (m_mode == TIE_FSI_FSI);

    // update Lagrange multipliers
    double normS1 = 0, normL1 = 0, normJ1 = 0;
    for (int np=0; np<2; ++np)
    {
        FETiedFluidFSISurface& s = (np == 0 ? m_ss : m_ms);
        for (int i=0; i<s.Elements(); ++i)
        {
            FESurfaceElement& el = s.Element(i);
            for (int j = 0; j<el.GaussPoints(); ++j)
            {
                FETiedFluidFSISurface::Data& d = static_cast<FETiedFluidFSISurface::Data&>(*el.GetMaterialPoint(j));

                if (d.m_pme) {
                    double epss = m_epss*d.m_epss;
                    double epst = m_epst*d.m_epst;
                    double epsn = m_epsn*d.m_epsn;

                    // Keep the reported tractions consistent with the multipliers about
                    // to be updated. NOTE: this has to happen BEFORE the updates below,
                    // otherwise the stored traction would be lambda_old + 2*eps*g. That
                    // stale value is what gets plotted on the augmentation that converges.
                    d.m_ts = d.m_Lmd + d.m_dg*epss;
                    // against a solid wall the traction is purely normal, exactly as in
                    // ProjectSurface; without the projection the value stored on the
                    // converging augmentation would keep a tangential residue of m_Lmt
                    d.m_tv = (bJ ? d.m_Lmt + d.m_gw*epst
                                 : d.m_nu*((d.m_Lmt*d.m_nu) + (d.m_gw*d.m_nu)*epst));
                    d.m_wn = (bJ ? d.m_Lmp + epsn*d.m_Jg : 0.0);

                    // solid traction multiplier, eq. (7.9.7)
                    d.m_Lmd = d.m_Lmd + d.m_dg*epss;
                    maxdg = max(maxdg, d.m_dg.norm());
                    normS1 += d.m_Lmd*d.m_Lmd;

                    // viscous traction multiplier. Against a solid wall the gap is
                    // g^w = -(w.n) n, so this drives the impermeability condition
                    // w.n = 0 and leaves the tangential slip unconstrained.
                    d.m_Lmt = d.m_Lmt + d.m_gw*epst;
                    maxgap = max(maxgap, d.m_gw.norm());
                    normL1 += d.m_Lmt*d.m_Lmt;

                    // normal velocity multiplier. There is no dilatation constraint
                    // against a solid wall, so there is no multiplier to update.
                    if (bJ) {
                        d.m_Lmp = d.m_Lmp + epsn*d.m_Jg;
                        maxJg = max(maxJg, fabs(d.m_Jg));
                        normJ1 += d.m_Lmp*d.m_Lmp;
                    }
                }
            }
        }
        if (m_btwo_pass == false) break;
    }
    normS1 = sqrt(normS1);
    normL1 = sqrt(normL1);
    normJ1 = sqrt(normJ1);

    // calculate relative norms
    double snorm = (normS1 != 0 ? fabs((normS1 - normS0) / normS1) : fabs(normS1 - normS0));
    double lnorm = (normL1 != 0 ? fabs((normL1 - normL0) / normL1) : fabs(normL1 - normL0));
    double pnorm = (normJ1 != 0 ? fabs((normJ1 - normJ0) / normJ1) : fabs(normJ1 - normJ0));

    // check convergence
    if ((m_gtol > 0) && (maxdg  > m_gtol)) bconv = false;
    if ((m_wtol > 0) && (maxgap > m_wtol)) bconv = false;
    if ((m_etol > 0) && (maxJg  > m_etol)) bconv = false;

    if ((m_atol > 0) && (snorm > m_atol)) bconv = false;
    if ((m_atol > 0) && (lnorm > m_atol)) bconv = false;
    if ((m_atol > 0) && (pnorm > m_atol)) bconv = false;

    if (naug < m_naugmin ) bconv = false;
    if (naug >= m_naugmax) bconv = true;

    feLog(" tied fluid-FSI interface # %d\n", GetID());
    feLog("                                CURRENT        REQUIRED\n");
    feLog("    solid multiplier      : %15le", snorm); if (m_atol > 0) feLog("%15le\n", m_atol); else feLog("       ***\n");
    feLog("    velocity multiplier   : %15le", lnorm); if (m_atol > 0) feLog("%15le\n", m_atol); else feLog("       ***\n");
    feLog("    dilatation multiplier : %15le", pnorm); if (m_atol > 0) feLog("%15le\n", m_atol); else feLog("       ***\n");
    feLog("    maximum gap           : %15le", maxdg);
    if (m_gtol > 0) feLog("%15le\n", m_gtol); else feLog("       ***\n");
    feLog("    max velocity gap%s: %15le", (m_mode == TIE_FSI_SOLID ? " (normal) " : "      "), maxgap);
    if (m_wtol > 0) feLog("%15le\n", m_wtol); else feLog("       ***\n");
    feLog("    maximum dilatation gap: %15le", maxJg);
    if (m_etol > 0) feLog("%15le\n", m_etol); else feLog("       ***\n");

    return bconv;
}

//-----------------------------------------------------------------------------
void FETiedFluidFSI::Serialize(DumpStream &ar)
{
    // store contact data
    FEContactInterface::Serialize(ar);

    // store contact surface data
    m_ss.Serialize(ar);
    m_ms.Serialize(ar);

	// serialize element pointers
	SerializeElementPointers(m_ss, m_ms, ar);
	SerializeElementPointers(m_ms, m_ss, ar);

    if (ar.IsShallow()) return;
    ar & m_dofU;
    ar & m_dofSU;
    ar & m_dofW;
    ar & m_dofEF;
    ar & m_mode;

    // A dump restart deserializes and resumes without re-running Init, so the owned
    // traction load has to be rebuilt here or the fluid traction would silently stop
    // being transferred to the solid after a restart. FESurface::Serialize re-finds
    // el.m_elem[] on loading, which is all Activate needs to rebuild its element list.
    if (ar.IsLoading())
    {
        delete m_pFSItrac;
        m_pFSItrac = nullptr;

        if (m_mode == TIE_FSI_SOLID)
        {
            m_pFSItrac = new FEFluidFSITraction(&ar.GetFEModel());
            m_pFSItrac->SetSurface(&m_ss);
            m_pFSItrac->SetShellBottom(m_bshellb);

            // NOTE: deliberately NOT calling Init() here. FESurfaceLoad::Init re-Inits
            // its surface, and FESurface::Init re-Inits every surface material point,
            // which would zero the gaps, the multipliers and m_pme that the two calls
            // above have just restored -- and Activate(), which is what would normally
            // re-project, does not run again on a restart. The tie would then transmit
            // nothing while the fluid traction kept being applied. Init() does nothing
            // else we need, so apply its two side effects directly instead.
            m_ss.SetInterfaceStatus(true);
            m_ss.SetShellBottom(m_bshellb);

            m_pFSItrac->Activate();
            m_pFSItrac->SetInterfaceOrientation();
        }
    }
}
