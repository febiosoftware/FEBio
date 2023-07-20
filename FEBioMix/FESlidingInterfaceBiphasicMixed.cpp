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
#include "FESlidingInterfaceBiphasicMixed.h"
#include "FEBiphasic.h"
#include "FECore/FEAnalysis.h"
#include "FECore/FENormalProjection.h"
#include <FECore/FELinearSystem.h>
#include "FECore/log.h"
#include <FEBioMech/FEBioMech.h>
#include "FEBioMix.h"
#include <FECore/FEModel.h>

//-----------------------------------------------------------------------------
// Define sliding interface parameters
BEGIN_FECORE_CLASS(FESlidingInterfaceBiphasicMixed, FEContactInterface)
	ADD_PARAMETER(m_atol     , "tolerance"          );
    ADD_PARAMETER(m_gtol     , "gaptol"             )->setUnits(UNIT_LENGTH);;
	ADD_PARAMETER(m_ptol     , "ptol"               );
	ADD_PARAMETER(m_epsn     , "penalty"            );
	ADD_PARAMETER(m_bautopen , "auto_penalty"       );
	ADD_PARAMETER(m_bupdtpen , "update_penalty"     );
	ADD_PARAMETER(m_btwo_pass, "two_pass"           );
	ADD_PARAMETER(m_knmult   , "knmult"             );
	ADD_PARAMETER(m_stol     , "search_tol"         );
	ADD_PARAMETER(m_epsp     , "pressure_penalty"   );
	ADD_PARAMETER(m_bsymm    , "symmetric_stiffness");
    ADD_PARAMETER(m_srad     , "search_radius"      )->setUnits(UNIT_LENGTH);;
	ADD_PARAMETER(m_nsegup   , "seg_up"             );
	ADD_PARAMETER(m_naugmin  , "minaug"             );
	ADD_PARAMETER(m_naugmax  , "maxaug"             );
	ADD_PARAMETER(m_breloc   , "node_reloc"         );
	ADD_PARAMETER(m_mu       , "fric_coeff"         );
	ADD_PARAMETER(m_phi      , "contact_frac"       );
	ADD_PARAMETER(m_bsmaug   , "smooth_aug"         );
    ADD_PARAMETER(m_bsmfls   , "smooth_fls"         );
    ADD_PARAMETER(m_bflips   , "flip_primary"       );
    ADD_PARAMETER(m_bflipm   , "flip_secondary"     );
    ADD_PARAMETER(m_bshellbs , "shell_bottom_primary"  );
    ADD_PARAMETER(m_bshellbm , "shell_bottom_secondary");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FESlidingSurfaceBiphasicMixed::Data::Data()
{
    m_Lmd   = 0.0;
    m_Lmt   = m_tr = vec3d(0,0,0);
    m_Lmp   = 0.0;
    m_epsn  = 1.0;
    m_epsp  = 1.0;
    m_pg    = 0.0;
    m_p1    = 0.0;
    m_mueff = 0.0;
    m_nu    = m_s1 = m_dg = vec3d(0,0,0);
    m_rs    = m_rsp = vec2d(0,0);
    m_bstick = false;
}

//-----------------------------------------------------------------------------
void FESlidingSurfaceBiphasicMixed::Data::Serialize(DumpStream& ar)
{
	FEBiphasicContactPoint::Serialize(ar);
	ar & m_dg;
    ar & m_Lmd;
    ar & m_Lmt;
    ar & m_epsn;
    ar & m_epsp;
    ar & m_p1;
    ar & m_nu;
	ar & m_s1;
    ar & m_tr;
	ar & m_rs;
	ar & m_rsp;
    ar & m_bstick;
}

//-----------------------------------------------------------------------------
// FESlidingSurfaceBiphasic
//-----------------------------------------------------------------------------

FESlidingSurfaceBiphasicMixed::FESlidingSurfaceBiphasicMixed(FEModel* pfem) : FEBiphasicContactSurface(pfem)
{
    m_bporo = false;
}

//-----------------------------------------------------------------------------
//! create material point data
FEMaterialPoint* FESlidingSurfaceBiphasicMixed::CreateMaterialPoint()
{
	return new FESlidingSurfaceBiphasicMixed::Data;
}

//-----------------------------------------------------------------------------
bool FESlidingSurfaceBiphasicMixed::Init()
{
	// get the displacement and fluid pressure variable indices.
	FEModel* fem = GetFEModel();
	m_varU = fem->GetDOFS().GetVariableIndex(FEBioMech::GetVariableName(FEBioMech::DISPLACEMENT)); assert(m_varU >= 0);
	m_varP = fem->GetDOFS().GetVariableIndex(FEBioMix::GetVariableName(FEBioMix::FLUID_PRESSURE)); assert(m_varP >= 0);

    // initialize surface data first
    if (FEBiphasicContactSurface::Init() == false) return false;
    
    // allocate node normals and contact tractions
    m_nn.assign(Nodes(), vec3d(0,0,0));
    m_tn.assign(Nodes(), vec3d(0,0,0));
    m_pn.assign(Nodes(), 0);
    
    // determine biphasic status
    m_poro.resize(Elements(),false);
    for (int i=0; i<Elements(); ++i)
    {
        // get the surface element
        FESurfaceElement& se = Element(i);
        
        // get the element this surface element belongs to
        FEElement* pe = se.m_elem[0];
        if (pe)
        {
            // get the material
            FEMaterial* pm = m_pfem->GetMaterial(pe->GetMatID());
            
            // see if this is a poro-elastic element
            FEBiphasic* biph = dynamic_cast<FEBiphasic*> (pm);
            if (biph) {
                m_poro[i] = true;
                m_bporo = true;
            }
        }
    }
    
    return true;
}

//-----------------------------------------------------------------------------
void FESlidingSurfaceBiphasicMixed::InitSlidingSurface()
{
    for (int i=0; i<Elements(); ++i)
    {
        FESurfaceElement& el = Element(i);
        int nint = el.GaussPoints();
        for (int j=0; j<nint; ++j)
        {
			Data& data = static_cast<Data&>(*el.GetMaterialPoint(j));
            // Store current surface projection values as previous
            data.m_rsp  = data.m_rs;
			data.m_pmep = data.m_pme;
        }
    }
}

//-----------------------------------------------------------------------------
//! Evaluate the nodal contact pressures by averaging values from surrounding
//! faces.  This function ensures that nodal contact pressures are always
//! positive, so that they can be used to detect free-draining status.

void FESlidingSurfaceBiphasicMixed::EvaluateNodalContactPressures()
{
    const int N = Nodes();
    
    // number of faces with non-zero contact pressure connected to this node
    vector<int> nfaces(N,0);
    
    // zero nodal contact pressures
    zero(m_pn);
    
    // loop over all elements
    for (int i=0; i<Elements(); ++i)
    {
        FESurfaceElement& el = Element(i);
        int ne = el.Nodes();
        
        // get the average contact pressure for that face
        double pn = 0;
        GetContactPressure(i, pn);
        
        if (pn > 0) {
            for (int j=0; j<ne; ++j)
            {
                m_pn[el.m_lnode[j]] += pn;
                ++nfaces[el.m_lnode[j]];
            }
        }
    }
    
    // get average over all contacting faces sharing that node
    for (int i=0; i<N; ++i)
        if (nfaces[i] > 0) m_pn[i] /= nfaces[i];
}

//-----------------------------------------------------------------------------
//! Evaluate the nodal contact tractions by averaging values from surrounding
//! faces.  This function ensures that nodal contact tractions are always
//! compressive, so that they can be used to detect free-draining status.

void FESlidingSurfaceBiphasicMixed::EvaluateNodalContactTractions()
{
    const int N = Nodes();
    
    // number of faces with non-zero contact pressure connected to this node
    vector<int> nfaces(N,0);
    
    // zero nodal contact tractions
    zero(m_tn);
    
    // loop over all elements
    for (int i=0; i<Elements(); ++i)
    {
        FESurfaceElement& el = Element(i);
        int ne = el.Nodes();
        
        // get the average contact traction and pressure for that face
        vec3d tn(0,0,0);
        GetContactTraction(i, tn);
        double pn = 0;
        GetContactPressure(i, pn);
        
        if (pn > 0) {
            for (int j=0; j<ne; ++j)
            {
                m_tn[el.m_lnode[j]] += tn;
                ++nfaces[el.m_lnode[j]];
            }
        }
    }
    
    // get average over all contacting faces sharing that node
    for (int i=0; i<N; ++i)
        if (nfaces[i] > 0) m_tn[i] /= nfaces[i];
}

//-----------------------------------------------------------------------------
//! This function calculates the node normal. Due to the piecewise continuity
//! of the surface elements this normal is not uniquely defined so in order to
//! obtain a unique normal the normal is averaged for each node over all the
//! element normals at the node

void FESlidingSurfaceBiphasicMixed::UpdateNodeNormals()
{
    const int MN = FEElement::MAX_NODES;
    vec3d y[MN];
    
    // zero nodal normals
    zero(m_nn);
    
    // loop over all elements
    for (int i=0; i<Elements(); ++i)
    {
        FESurfaceElement& el = Element(i);
        int ne = el.Nodes();
        
        // get the nodal coordinates
        for (int j=0; j<ne; ++j) y[j] = Node(el.m_lnode[j]).m_rt;
        
        // calculate the normals
        for (int j=0; j<ne; ++j)
        {
            int jp1 = (j+1)%ne;
            int jm1 = (j+ne-1)%ne;
            vec3d n = (y[jp1] - y[j]) ^ (y[jm1] - y[j]);
            m_nn[el.m_lnode[j]] += n;
        }
    }
    
    // normalize all vectors
    const int N = Nodes();
    for (int i=0; i<N; ++i) m_nn[i].unit();
}

//-----------------------------------------------------------------------------
vec3d FESlidingSurfaceBiphasicMixed::GetContactForce()
{
    return m_Ft;
}

//-----------------------------------------------------------------------------
double FESlidingSurfaceBiphasicMixed::GetContactArea()
{
    // initialize contact area
    double a = 0;
    
    // loop over all elements of the primary surface
    for (int n=0; n<Elements(); ++n)
    {
        FESurfaceElement& el = Element(n);
        int nint = el.GaussPoints();
        
        // evaluate the contact force for that element
        for (int i=0; i<nint; ++i)
        {
            // get data for this integration point
            Data& data = static_cast<Data&>(*el.GetMaterialPoint(i));
            if (data.m_Ln > 0)
            {
                // get the base vectors
                vec3d g[2];
                CoBaseVectors(el, i, g);
                
                // normal (magnitude = area)
                vec3d n = g[0] ^ g[1];
                
                // gauss weight
                double w = el.GaussWeights()[i];
                
                // contact force
                a += n.norm()*w;
            }
        }
    }
    
    return a;
}

//-----------------------------------------------------------------------------
vec3d FESlidingSurfaceBiphasicMixed::GetFluidForce()
{
	DOFS& dofs = GetFEModel()->GetDOFS();
	int degree_d = dofs.GetVariableInterpolationOrder(m_varU);
	int degree_p = dofs.GetVariableInterpolationOrder(m_varP);

    const int MN = FEElement::MAX_NODES;
    double pn[MN];
    
    // initialize contact force
    vec3d f(0,0,0);
    
    // loop over all elements of the surface
    for (int n=0; n<Elements(); ++n)
    {
        FESurfaceElement& el = Element(n);
        int nseln = el.Nodes();
        
        // nodal pressures
		int npdof = el.ShapeFunctions(degree_p);
        for (int i=0; i<npdof; ++i) pn[i] = GetMesh()->Node(el.m_node[i]).get(m_dofP);
        
       
        // evaluate the fluid force for that element
		int nint = el.GaussPoints();
		for (int i=0; i<nint; ++i)
        {
            // get the base vectors
            vec3d g[2];
            CoBaseVectors(el, i, g);
            // normal (magnitude = area)
            vec3d n = g[0] ^ g[1];
            // gauss weight
            double w = el.GaussWeights()[i];
            // fluid pressure
            double p = el.eval(degree_p, pn, i);
            // contact force
            f += n*(w*p);
        }
    }
    
    return f;
}

//-----------------------------------------------------------------------------
void FESlidingSurfaceBiphasicMixed::Serialize(DumpStream& ar)
{
	FEBiphasicContactSurface::Serialize(ar);
	ar & m_bporo;
	ar & m_poro;
	ar & m_nn;
	ar & m_pn;
	ar & m_tn;
	ar & m_Ft;
}

//-----------------------------------------------------------------------------
void FESlidingSurfaceBiphasicMixed::GetVectorGap(int nface, vec3d& pg)
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
void FESlidingSurfaceBiphasicMixed::GetContactPressure(int nface, double& pg)
{
    FESurfaceElement& el = Element(nface);
    int ni = el.GaussPoints();
    pg = 0;
	for (int k = 0; k < ni; ++k)
	{
		Data& data = static_cast<Data&>(*el.GetMaterialPoint(k));
		pg += data.m_Ln;
	}
    pg /= ni;
}

//-----------------------------------------------------------------------------
void FESlidingSurfaceBiphasicMixed::GetContactTraction(int nface, vec3d& pt)
{
    FESurfaceElement& el = Element(nface);
    int ni = el.GaussPoints();
    pt = vec3d(0,0,0);
	for (int k = 0; k < ni; ++k)
	{
		Data& data = static_cast<Data&>(*el.GetMaterialPoint(k));
		pt += data.m_tr;
	}
    pt /= ni;
}

//-----------------------------------------------------------------------------
void FESlidingSurfaceBiphasicMixed::GetSlipTangent(int nface, vec3d& pt)
{
    FESurfaceElement& el = Element(nface);
    int ni = el.GaussPoints();
    pt = vec3d(0,0,0);
	for (int k = 0; k < ni; ++k)
	{
		Data& data = static_cast<Data&>(*el.GetMaterialPoint(k));
		if (!data.m_bstick) pt += data.m_s1;
	}
    pt /= ni;
}

//-----------------------------------------------------------------------------
void FESlidingSurfaceBiphasicMixed::GetMuEffective(int nface, double& pg)
{
    FESurfaceElement& el = Element(nface);
    int ni = el.GaussPoints();
    pg = 0;
	for (int k = 0; k < ni; ++k)
	{
		Data& data = static_cast<Data&>(*el.GetMaterialPoint(k));
		pg += data.m_mueff;
	}
    pg /= ni;
}

//-----------------------------------------------------------------------------
void FESlidingSurfaceBiphasicMixed::GetLocalFLS(int nface, double& pg)
{
    FESurfaceElement& el = Element(nface);
    int ni = el.GaussPoints();
    pg = 0;
    for (int k = 0; k < ni; ++k)
    {
        Data& data = static_cast<Data&>(*el.GetMaterialPoint(k));
        pg += data.m_fls;
    }
    pg /= ni;
}

//-----------------------------------------------------------------------------
void FESlidingSurfaceBiphasicMixed::GetNodalVectorGap(int nface, vec3d* pg)
{
    FESurfaceElement& el = Element(nface);
    int ni = el.GaussPoints();
    vec3d gi[FEElement::MAX_INTPOINTS];
	for (int k = 0; k < ni; ++k)
	{
		Data& data = static_cast<Data&>(*el.GetMaterialPoint(k));
		gi[k] = data.m_dg;
	}
    el.project_to_nodes(gi, pg);
}

//-----------------------------------------------------------------------------
void FESlidingSurfaceBiphasicMixed::GetNodalContactPressure(int nface, double* pg)
{
	DOFS& dofs = GetFEModel()->GetDOFS();
	int degree_p = dofs.GetVariableInterpolationOrder(m_varP);

    FESurfaceElement& el = Element(nface);
	int npdof = el.ShapeFunctions(degree_p);
    for (int k=0; k<npdof; ++k)
        pg[k] = m_pn[el.m_lnode[k]];
}

//-----------------------------------------------------------------------------
void FESlidingSurfaceBiphasicMixed::GetStickStatus(int nface, double& pg)
{
    FESurfaceElement& el = Element(nface);
    int ni = el.GaussPoints();
    pg = 0;
	for (int k = 0; k < ni; ++k)
	{
		Data& data = static_cast<Data&>(*el.GetMaterialPoint(k));
		if (data.m_bstick) pg += 1.0;
	}
    pg /= ni;
}

//-----------------------------------------------------------------------------
void FESlidingSurfaceBiphasicMixed::GetNodalContactTraction(int nface, vec3d* tn)
{
    FESurfaceElement& el = Element(nface);
    for (int k=0; k<el.Nodes(); ++k)
        tn[k] = m_tn[el.m_lnode[k]];
}

//-----------------------------------------------------------------------------
void FESlidingSurfaceBiphasicMixed::UnpackLM(FEElement& el, vector<int>& lm)
{
	DOFS& dofs = GetFEModel()->GetDOFS();
	int degree_d = dofs.GetVariableInterpolationOrder(m_varU);
	int degree_p = dofs.GetVariableInterpolationOrder(m_varP);

	// We should allocate the number of shapefunctions for each variable,
	// but for now, we allocate assuming the shape functions equals the nodes.
	// We can do this, because the nodes that don't have a pressure dof,
	// will have their id set to -1. 
	int N = el.Nodes();
	lm.assign(N * 4, -1);

	// pack the equation numbers
	for (int i = 0; i<N; ++i)
	{
		int n = el.m_node[i];

		FENode& node = m_pMesh->Node(n);
		vector<int>& id = node.m_ID;

		// first the displacement dofs
		lm[3 * i    ] = id[m_dofX];
		lm[3 * i + 1] = id[m_dofY];
		lm[3 * i + 2] = id[m_dofZ];

		// now the pressure dofs
		if (m_dofP >= 0) lm[3 * N + i] = id[m_dofP];
	}
}

//-----------------------------------------------------------------------------
// FESlidingInterfaceBiphasicMixed
//-----------------------------------------------------------------------------

FESlidingInterfaceBiphasicMixed::FESlidingInterfaceBiphasicMixed(FEModel* pfem) : FEContactInterface(pfem), m_ss(pfem), m_ms(pfem)
{
    static int count = 1;
    SetID(count++);
    
    // initial values
    m_knmult = 0;
    m_atol = 0.1;
    m_epsn = 1;
    m_epsp = 1;
    m_btwo_pass = false;
    m_stol = 0.01;
    m_bsymm = true;
    m_srad = 1.0;
    m_gtol = 0;
    m_ptol = 0;
    m_nsegup = 0;
    m_bautopen = false;
    m_breloc = false;
    m_bsmaug = false;
    m_bsmfls = true;
    m_bupdtpen = false;
    m_mu = 0.0;
    m_phi = 0.0;
    
    m_naugmin = 0;
    m_naugmax = 10;
    
    m_bfreeze = false;
    m_bflipm = m_bflips = false;
    m_bshellbm = m_bshellbs = false;

    m_dofP = pfem->GetDOFIndex("p");
    
    // set parents
    m_ss.SetContactInterface(this);
    m_ms.SetContactInterface(this);

    m_ss.SetSibling(&m_ms);
    m_ms.SetSibling(&m_ss);
}

//-----------------------------------------------------------------------------

FESlidingInterfaceBiphasicMixed::~FESlidingInterfaceBiphasicMixed()
{
}

//-----------------------------------------------------------------------------
bool FESlidingInterfaceBiphasicMixed::Init()
{
    // initialize surface data
    if (m_ss.Init() == false) return false;
    if (m_ms.Init() == false) return false;
    
    // Flip secondary and primary surfaces, if requested.
    // Note that we turn off those flags because otherwise we keep flipping, each time we get here (e.g. in optimization)
    // TODO: Of course, we shouldn't get here more than once. I think we also get through the FEModel::Reset, so I'll have
    //       look into that.
    if (m_bflips) { m_ss.Invert(); m_bflips = false; }
    if (m_bflipm) { m_ms.Invert(); m_bflipm = false; }
    if (m_bshellbs) { m_ss.SetShellBottom(m_bshellbs); m_bshellbs = false; }
    if (m_bshellbm) { m_ms.SetShellBottom(m_bshellbm); m_bshellbm = false; }
    
    return true;
}

//-----------------------------------------------------------------------------
void FESlidingInterfaceBiphasicMixed::BuildMatrixProfile(FEGlobalMatrix& K)
{
    FEModel& fem = *GetFEModel();
    FEMesh& mesh = fem.GetMesh();
    
    // get the DOFS
    const int dof_X = fem.GetDOFIndex("x");
    const int dof_Y = fem.GetDOFIndex("y");
    const int dof_Z = fem.GetDOFIndex("z");
    const int dof_P = fem.GetDOFIndex("p");
    const int dof_RU = fem.GetDOFIndex("Ru");
    const int dof_RV = fem.GetDOFIndex("Rv");
    const int dof_RW = fem.GetDOFIndex("Rw");
    
    vector<int> lm(7*FEElement::MAX_NODES*2);
    
    int npass = (m_btwo_pass?2:1);
    for (int np=0; np<npass; ++np)
    {
        FESlidingSurfaceBiphasicMixed& ss = (np == 0? m_ss : m_ms);
        
        int k, l;
        for (int j=0; j<ss.Elements(); ++j)
        {
            FESurfaceElement& se = ss.Element(j);
            int nint = se.GaussPoints();
            int* sn = &se.m_node[0];
            for (k=0; k<nint; ++k)
            {
				FESlidingSurfaceBiphasicMixed::Data& pt = static_cast<FESlidingSurfaceBiphasicMixed::Data&>(*se.GetMaterialPoint(k));
                FESurfaceElement* pe = pt.m_pme;
                if (pe != 0)
                {
                    FESurfaceElement& me = *pe;
                    int* mn = &me.m_node[0];
                    
                    assign(lm, -1);
                    
                    int nseln = se.Nodes();
                    int nmeln = me.Nodes();
                    
                    for (l=0; l<nseln; ++l)
                    {
                        vector<int>& id = mesh.Node(sn[l]).m_ID;
                        lm[7*l  ] = id[dof_X];
                        lm[7*l+1] = id[dof_Y];
                        lm[7*l+2] = id[dof_Z];
                        lm[7*l+3] = id[dof_P];
                        lm[7*l+4] = id[dof_RU];
                        lm[7*l+5] = id[dof_RV];
                        lm[7*l+6] = id[dof_RW];
                    }
                    
                    for (l=0; l<nmeln; ++l)
                    {
                        vector<int>& id = mesh.Node(mn[l]).m_ID;
                        lm[7*(l+nseln)  ] = id[dof_X];
                        lm[7*(l+nseln)+1] = id[dof_Y];
                        lm[7*(l+nseln)+2] = id[dof_Z];
                        lm[7*(l+nseln)+3] = id[dof_P];
                        lm[7*(l+nseln)+4] = id[dof_RU];
                        lm[7*(l+nseln)+5] = id[dof_RV];
                        lm[7*(l+nseln)+6] = id[dof_RW];
                    }
                    
                    K.build_add(lm);
                }
            }
        }
    }
}

//-----------------------------------------------------------------------------
void FESlidingInterfaceBiphasicMixed::UpdateAutoPenalty()
{
    // calculate the penalty
    if (m_bautopen)
    {
        CalcAutoPenalty(m_ss);
        CalcAutoPenalty(m_ms);
        CalcAutoPressurePenalty(m_ss);
        CalcAutoPressurePenalty(m_ms);
    }
}

//-----------------------------------------------------------------------------
//! This function is called during the initialization
void FESlidingInterfaceBiphasicMixed::Activate()
{
    // don't forget to call base member
    FEContactInterface::Activate();
    
    UpdateAutoPenalty();
    
    // update sliding interface data
    Update();
}

//-----------------------------------------------------------------------------
void FESlidingInterfaceBiphasicMixed::CalcAutoPenalty(FESlidingSurfaceBiphasicMixed& s)
{
    // loop over all surface elements
    for (int i=0; i<s.Elements(); ++i)
    {
        // get the surface element
        FESurfaceElement& el = s.Element(i);
        
        // calculate a penalty
        double eps = AutoPenalty(el, s);
        
        // assign to integation points of surface element
        int nint = el.GaussPoints();
        for (int j=0; j<nint; ++j)
        {
			FESlidingSurfaceBiphasicMixed::Data& pt = static_cast<FESlidingSurfaceBiphasicMixed::Data&>(*el.GetMaterialPoint(j));
			pt.m_epsn = eps;
        }
    }
}

//-----------------------------------------------------------------------------
void FESlidingInterfaceBiphasicMixed::CalcAutoPressurePenalty(FESlidingSurfaceBiphasicMixed& s)
{
    // loop over all surface elements
    for (int i=0; i<s.Elements(); ++i)
    {
        // get the surface element
        FESurfaceElement& el = s.Element(i);
        
        // calculate a penalty
        double eps = AutoPressurePenalty(el, s);
        
        // assign to integation points of surface element
        int nint = el.GaussPoints();
        for (int j=0; j<nint; ++j)
        {
			FESlidingSurfaceBiphasicMixed::Data& pt = static_cast<FESlidingSurfaceBiphasicMixed::Data&>(*el.GetMaterialPoint(j));
			pt.m_epsp = eps;
        }
    }
}

//-----------------------------------------------------------------------------

double FESlidingInterfaceBiphasicMixed::AutoPressurePenalty(FESurfaceElement& el, FESlidingSurfaceBiphasicMixed& s)
{
    // get the mesh
    FEMesh& m = GetFEModel()->GetMesh();
    
    // evaluate element surface normal at parametric center
    vec3d t[2];
    s.CoBaseVectors0(el, 0, 0, t);
    vec3d n = t[0] ^ t[1];
    n.unit();
    
    // get the element this surface element belongs to
    FEElement* pe = el.m_elem[0];
    if (pe == 0) return 0.0;

    // get the material
    FEMaterial* pm = GetFEModel()->GetMaterial(pe->GetMatID());
        
    // see if this is a poro-elastic element
    FEBiphasic* biph = dynamic_cast<FEBiphasic*> (pm);
    if (biph == 0) return 0.0;

    // get a material point
    FEMaterialPoint& mp = *pe->GetMaterialPoint(0);
    FEElasticMaterialPoint& ept = *(mp.ExtractData<FEElasticMaterialPoint>());
            
    // setup the material point
    ept.m_F = mat3dd(1.0);
    ept.m_J = 1;
    ept.m_s.zero();
            
    // if this is a poroelastic element, then get the permeability tensor
    FEBiphasicMaterialPoint& pt = *(mp.ExtractData<FEBiphasicMaterialPoint>());
    pt.m_p = 0;
    pt.m_w = vec3d(0,0,0);
            
    double K[3][3];
    biph->Permeability(K, mp);
            
    double eps = n.x*(K[0][0]*n.x+K[0][1]*n.y+K[0][2]*n.z)
    +n.y*(K[1][0]*n.x+K[1][1]*n.y+K[1][2]*n.z)
    +n.z*(K[2][0]*n.x+K[2][1]*n.y+K[2][2]*n.z);
    
	// get the area of the surface element
	double A = s.FaceArea(el);

	// get the volume of the volume element
	double V = m.ElementVolume(*pe);

    return eps*A/V;
}

//-----------------------------------------------------------------------------
void FESlidingInterfaceBiphasicMixed::ProjectSurface(FESlidingSurfaceBiphasicMixed& ss, FESlidingSurfaceBiphasicMixed& ms, bool bupseg, bool bmove)
{
	DOFS& dofs = GetFEModel()->GetDOFS();
	int degree_d = dofs.GetVariableInterpolationOrder(ss.m_varU);
	int degree_p = dofs.GetVariableInterpolationOrder(ss.m_varP);

    FEMesh& mesh = GetFEModel()->GetMesh();
    
    // initialize projection data
    FENormalProjection np(ms);
    np.SetTolerance(m_stol);
    np.SetSearchRadius(m_srad);
    np.Init();
    double psf = GetPenaltyScaleFactor();

    // if we need to project the nodes onto the secondary surface,
    // let's do this first
    if (bmove)
    {
        int NN = ss.Nodes();
        int NE = ss.Elements();
        // first we need to calculate the node normals
        vector<vec3d> normal; normal.assign(NN, vec3d(0,0,0));
        for (int i=0; i<NE; ++i)
        {
            FESurfaceElement& el = ss.Element(i);
            int ne = el.Nodes();
            for (int j=0; j<ne; ++j)
            {
                vec3d r0 = ss.Node(el.m_lnode[ j         ]).m_rt;
                vec3d rp = ss.Node(el.m_lnode[(j+   1)%ne]).m_rt;
                vec3d rm = ss.Node(el.m_lnode[(j+ne-1)%ne]).m_rt;
                vec3d n = (rp - r0)^(rm - r0);
                normal[el.m_lnode[j]] += n;
            }
        }
        for (int i=0; i<NN; ++i) normal[i].unit();
        
        // loop over all nodes
        for (int i=0; i<NN; ++i)
        {
            FENode& node = ss.Node(i);
            
            // get the spatial nodal coordinates
            vec3d rt = node.m_rt;
            vec3d nu = normal[i];
            
            // project onto the secondary surface
            vec3d q;
            double rs[2] = {0,0};
            FESurfaceElement* pme = np.Project(rt, nu, rs);
            if (pme)
            {
                // the node could potentially be in contact
                // find the global location of the intersection point
                vec3d q = ms.Local2Global(*pme, rs[0], rs[1]);
                
                // calculate the gap function
                // NOTE: this has the opposite sign compared
                // to Gerard's notes.
                double gap = nu*(rt - q);
                
                if (gap>0) node.m_r0 = node.m_rt = q;
            }
        }
    }
    
    // loop over all integration points
#pragma omp parallel for
    for (int i=0; i<ss.Elements(); ++i)
    {
        FESurfaceElement& el = ss.Element(i);
        bool sporo = ss.m_poro[i];
        
        int ne = el.Nodes();
        int nint = el.GaussPoints();
        double ps[FEElement::MAX_INTPOINTS];
        // get the nodal pressures
        if (sporo)
        {
			int npd = el.ShapeFunctions(degree_p);
            for (int j=0; j<npd; ++j) ps[j] = mesh.Node(el.m_node[j]).get(m_dofP);
        }
        
        for (int j=0; j<nint; ++j)
        {
            // get the integration point data
			FESlidingSurfaceBiphasicMixed::Data& pt = static_cast<FESlidingSurfaceBiphasicMixed::Data&>(*el.GetMaterialPoint(j));

            // calculate the global position of the integration point
            vec3d r = ss.Local2Global(el, j);
            
            // get the pressure at the integration point
            double p1 = 0;
            if (sporo) p1 = el.eval(degree_p, ps, j);

            // calculate the normal at this integration point
            vec3d nu = ss.SurfaceNormal(el, j);
            
            // first see if the old intersected face is still good enough
            FESurfaceElement* pme = pt.m_pme;
            double rs[2] = {0,0};
            if (pme)
            {
                double g;
                
                // see if the ray intersects this element
                if (ms.Intersect(*pme, r, nu, rs, g, m_stol))
                {
                    pt.m_rs[0] = rs[0];
                    pt.m_rs[1] = rs[1];
                }
                else
                {
                    pme = 0;
                }
            }
            
            // find the intersection point with the secondary surface
            if (pme == 0 && bupseg) pme = np.Project(r, nu, rs);
            
            pt.m_pme = pme;
            pt.m_nu = nu;
            pt.m_rs[0] = rs[0];
            pt.m_rs[1] = rs[1];
            if (pme)
            {
                // the node could potentially be in contact
                // find the global location of the intersection point
                vec3d q = ms.Local2Global(*pme, rs[0], rs[1]);
                
                // calculate the gap function
                // NOTE: this has the opposite sign compared
                // to Gerard's notes.
                double g = nu*(r - q);
                
                double eps = m_epsn*pt.m_epsn*psf;
                
                double Ln = pt.m_Lmd + eps*g;
                
                pt.m_gap = (g <= m_srad? g : 0);
                
                // calculate the pressure gap function
                bool mporo = ms.m_poro[pme->m_lid];
                
                if ((Ln >= 0) && (g <= m_srad))
                {
                    
                    // get the pressure at the projection point
                    double p2 = 0;
                    if (mporo) {
                        double pm[FEElement::MAX_NODES];
						int npdof = pme->ShapeFunctions(degree_p);
                        for (int k=0; k<npdof; ++k) pm[k] = mesh.Node(pme->m_node[k]).get(m_dofP);
                        p2 = pme->eval(degree_p, pm, rs[0], rs[1]);
                    }
                    if (sporo) {
                        pt.m_p1 = p1;
                        if (mporo) {
                            pt.m_pg = p1 - p2;
                        }
                    }
                    else if (mporo) {
                        pt.m_p1 = p2;
                    }
                }
                else
                {
                    pt.m_Lmd = 0;
                    pt.m_pme = 0;
                    pt.m_gap = 0;
                    pt.m_dg = pt.m_Lmt = vec3d(0,0,0);
                    if (sporo || mporo) {
                        pt.m_Lmp = 0;
                        pt.m_pg = 0;
                        pt.m_p1 = 0;
                    }
                }
            }
            else
            {
                // the node is not in contact
                pt.m_Lmd = 0;
                pt.m_gap = 0;
                pt.m_dg = pt.m_Lmt = vec3d(0,0,0);
                if (sporo) {
                    pt.m_Lmp = 0;
                    pt.m_pg = 0;
                    pt.m_p1 = 0;
                }
            }
        }
    }
}

//-----------------------------------------------------------------------------

void FESlidingInterfaceBiphasicMixed::Update()
{
	DOFS& dofs = GetFEModel()->GetDOFS();
	int degree_p = dofs.GetVariableInterpolationOrder(m_ss.m_varP);

    static int naug = 0;
    static int biter = 0;
    
    FEModel& fem = *GetFEModel();
    
    // get the iteration number
    // we need this number to see if we can do segment updates or not
    // also reset number of iterations after each augmentation
    FEAnalysis* pstep = fem.GetCurrentStep();
    FESolver* psolver = pstep->GetFESolver();
    if (psolver->m_niter == 0) {
        biter = 0;
        naug = psolver->m_naug;
        // check update of auto-penalty
        if (m_bupdtpen) UpdateAutoPenalty();
    } else if (psolver->m_naug > naug) {
        biter = psolver->m_niter;
        naug = psolver->m_naug;
    }
    int niter = psolver->m_niter - biter;
    bool bupseg = ((m_nsegup == 0)? true : (niter <= m_nsegup));
    // get the logfile
    //	Logfile& log = GetLogfile();
    //	log.printf("seg_up iteration # %d\n", niter+1);
    
    // project the surfaces onto each other
    // this will update the gap functions as well
    static bool bfirst = true;
    ProjectSurface(m_ss, m_ms, bupseg, (m_breloc && bfirst));
    if (m_btwo_pass || m_ms.m_bporo) ProjectSurface(m_ms, m_ss, bupseg);
    bfirst = false;
    
    // Call InitSlidingSurface on the first iteration of each time step
	int nsolve_iter = psolver->m_niter;
    if (nsolve_iter == 0)
    {
        m_ss.InitSlidingSurface();
        if (m_btwo_pass) m_ms.InitSlidingSurface();
        m_bfreeze = false;
    }
    
    // Update the net contact pressures
    UpdateContactPressures();
    
    if (niter == 0) m_bfreeze = false;
    
    // set poro flag
    bool bporo = (m_ss.m_bporo || m_ms.m_bporo);
    
    // only continue if we are doing a poro-elastic simulation
    if (bporo == false) return;
    
    // update node normals
    m_ss.UpdateNodeNormals();
    if (bporo) m_ms.UpdateNodeNormals();
    
    // Now that the nodes have been projected, we need to figure out
    // if we need to modify the constraints on the pressure dofs.
    // If the nodes are not in contact, they must be free
    // draining. Since all nodes have been previously marked to be
    // free-draining in MarkFreeDraining(), we just need to reverse
    // this setting here, for nodes that are in contact.
    
    // Next, we loop over each surface, visiting the nodes
    // and finding out if that node is in contact or not
    int npass = (m_btwo_pass?2:1);
    for (int np=0; np<npass; ++np)
    {
		FESlidingSurfaceBiphasicMixed& ss = (np == 0? m_ss : m_ms);
		FESlidingSurfaceBiphasicMixed& ms = (np == 0? m_ms : m_ss);
        
        // loop over all the nodes of the primary surface
        for (int n=0; n<ss.Nodes(); ++n) {
            FENode& node = ss.Node(n);
            int id = node.m_ID[m_dofP];
            if ((id < -1) && (ss.m_pn[n] > 0))
            {
                // mark node as non-free-draining (= pos ID)
                node.m_ID[m_dofP] = -id-2;
            }
        }
        
        // loop over all nodes of the secondary surface
        // the secondary surface is trickier since we need
        // to look at the primary surface's projection
        if (ms.m_bporo) {
            FENormalProjection np(ss);
            np.SetTolerance(m_stol);
            np.SetSearchRadius(m_srad);
            np.Init();
            
            for (int n=0; n<ms.Nodes(); ++n)
            {
                // get the node
                FENode& node = ms.Node(n);
                
                // project it onto the primary surface
                double rs[2] = {0,0};
                FESurfaceElement* pse = np.Project(node.m_rt, ms.m_nn[n], rs);
                
                if (pse)
                {
                    // we found an element, so let's see if it's even remotely close to contact
                    // find the global location of the intersection point
                    vec3d q = ss.Local2Global(*pse, rs[0], rs[1]);
                    
                    // calculate the gap function
                    double g = ms.m_nn[n]*(node.m_rt - q);
                    
                    if (fabs(g) <= m_srad)
                    {
                        // we found an element so let's calculate the nodal traction values for this element
                        // get the normal tractions at the nodes
                        
                        double tn[FEElement::MAX_NODES];
						int N = pse->ShapeFunctions(degree_p);
                        for (int i=0; i<N; ++i)
                            tn[i] = ss.m_pn[pse->m_lnode[i]];
                        
                        // now evaluate the traction at the intersection point
                        double tp = pse->eval(degree_p, tn, rs[0], rs[1]);
                        
                        // if tp > 0, mark node as non-free-draining. (= pos ID)
                        int id = node.m_ID[m_dofP];
                        if ((id < -1) && (tp > 0))
                        {
                            // mark as non free-draining
                            node.m_ID[m_dofP] = -id-2;
                        }
                    }
                }
            }
        }
    }
}

//-----------------------------------------------------------------------------
vec3d FESlidingInterfaceBiphasicMixed::SlipTangent(FESlidingSurfaceBiphasicMixed& ss, const int nel, const int nint, FESlidingSurfaceBiphasicMixed& ms, double& dh, vec3d& r)
{
    vec3d s1(0,0,0);
    dh = 0;
    r = vec3d(0,0,0);
    
    // get primary surface element
    FESurfaceElement& se = ss.Element(nel);
    
    // get integration point data
	FESlidingSurfaceBiphasicMixed::Data& data = static_cast<FESlidingSurfaceBiphasicMixed::Data&>(*se.GetMaterialPoint(nint));
	double g = data.m_gap;
    vec3d nu = data.m_nu;
    
    // find secondary surface element
    FESurfaceElement* pme = data.m_pme;
    
    // calculate previous positions
    vec3d x2p = ms.Local2GlobalP(*pme, data.m_rs[0], data.m_rs[1]);
    vec3d x1p = ss.Local2GlobalP(se, nint);
    
    // calculate dx2
    vec3d x2 = ms.Local2Global(*pme, data.m_rs[0], data.m_rs[1]);
    vec3d dx2 = x2 - x2p;
    
    // calculate dx1
    vec3d x1 = ss.Local2Global(se, nint);
    vec3d dx1 = x1 - x1p;
    
    // get current and previous covariant basis vectors
    vec3d gscov[2], gscovp[2];
    ss.CoBaseVectors(se, nint, gscov);
    ss.CoBaseVectorsP(se, nint, gscovp);
    
    // calculate delta gscov
    vec3d dgscov[2];
    dgscov[0] = gscov[0] - gscovp[0];
    dgscov[1] = gscov[1] - gscovp[1];
    
    // calculate m, J, Nhat
    vec3d m = ((dgscov[0] ^ gscov[1]) + (gscov[0] ^ dgscov[1]));
    double detJ = (gscov[0] ^ gscov[1]).norm();
    mat3d Nhat = (mat3dd(1) - (nu & nu));
    
    // calculate c
    vec3d c = Nhat*m*(1.0/detJ);
    
    // calculate slip direction s1
    double norm = (Nhat*(c*(-g) + dx1 - dx2)).norm();
    if (norm != 0)
    {
        s1 = (Nhat*(c*(-g) + dx1 - dx2))/norm;
        dh = norm;
        r = c*(-g) + dx1 - dx2;
    }
    
    return s1;
    
}

//-----------------------------------------------------------------------------
vec3d FESlidingInterfaceBiphasicMixed::ContactTraction(FESlidingSurfaceBiphasicMixed& ss, const int nel, const int n, FESlidingSurfaceBiphasicMixed& ms, double& pn)
{
	DOFS& dofs = GetFEModel()->GetDOFS();
	int degree_d = dofs.GetVariableInterpolationOrder(ss.m_varU);
	int degree_p = dofs.GetVariableInterpolationOrder(ss.m_varP);

    vec3d s1(0,0,0);
    vec3d dr(0,0,0);
    vec3d t(0,0,0);
    pn = 0;
    double tn = 0, ts = 0, mueff = 0;
    double psf = GetPenaltyScaleFactor();

    // get the mesh
    FEMesh& m = GetFEModel()->GetMesh();

	// get the primary surface element
	FESurfaceElement& se = ss.Element(nel);

    // get the integration point data
	FESlidingSurfaceBiphasicMixed::Data& data = static_cast<FESlidingSurfaceBiphasicMixed::Data&>(*se.GetMaterialPoint(n));

    // penalty
    double eps = m_epsn*data.m_epsn*psf;
    
    // normal gap
    double g = data.m_gap;
    
    // normal traction Lagrange multiplier
    double Lm = data.m_Lmd;
    
    // vector traction Lagrange multiplier
    vec3d Lt = data.m_Lmt;
    
    // get the normal at this integration point
    vec3d nu = data.m_nu;
    
    // get the fluid pressure at this integration point
    double p = data.m_p1;
    
    // get poro status of primary surface
    bool sporo = ss.m_poro[nel];
    
    // get current and previous secondary elements
    FESurfaceElement* pme = data.m_pme;
    FESurfaceElement* pmep = data.m_pmep;
    
    // zero the effective friction coefficient
    data.m_mueff = 0.0;
    data.m_fls = 0.0;
    data.m_s1 = vec3d(0,0,0);
    
    // get local FLS from element projection
    double fls = 0;
    if (m_bsmfls) {
        double lfls[FEElement::MAX_INTPOINTS];
        ss.GetGPLocalFLS(nel, lfls);
        fls = lfls[n];
    }
    
    // if we just returned from an augmentation, do not update stick or slip status
    if (m_bfreeze && pme) {
        if (data.m_bstick) {
            // calculate current global position of the integration point
            vec3d xo = ss.Local2Global(se, n);
            
            // calculate current global position of the previous intersection point
            vec3d xt = ms.Local2Global(*pmep, data.m_rsp[0], data.m_rsp[1]);
            
            // calculate vector gap
            vec3d dg = xt - xo;
            
            // calculate trial stick traction, normal component, shear component
            t = Lt + dg*eps;
            tn = t*nu;
            ts = (t - nu*tn).norm();
            
            // contact pressure
            pn = MBRACKET(-tn);
            
            // calculate effective friction coefficient
            if (pn > 0)
            {
                data.m_mueff = ts/pn;
                data.m_fls = m_bsmfls ? fls : p/pn;
            }

            // store the previous values as the current
            data.m_pme = data.m_pmep;
            data.m_rs = data.m_rsp;
            
            // recalculate gap
            data.m_dg = dg;
            
            // recalculate pressure gap
            bool mporo = ms.m_poro[pme->m_lid];
            if (sporo && mporo)
            {
                double pm[FEElement::MAX_NODES];
				int npdof = pme->ShapeFunctions(degree_p);
                for (int k=0; k<npdof; ++k) pm[k] = m.Node(pme->m_node[k]).get(m_dofP);
                double p2 = pme->eval(degree_p, pm, data.m_rs[0], data.m_rs[1]);
                data.m_pg = p - p2;
            }
            
        }
        else {
            // recalculate contact pressure for slip
            pn = MBRACKET(Lm + eps*g);
            
            if (pn != 0)
            {
                
                double dh = 0;
                
                // slip direction
                s1 = SlipTangent(ss, nel, n, ms, dh, dr);
                
                // calculate effective friction coefficient
                data.m_fls = m_bsmfls ? fls : p/pn;
                data.m_mueff = m_mu*(1.0-(1.0-m_phi)*data.m_fls);
                data.m_mueff = MBRACKET(data.m_mueff);

                // total traction
                t = nu*(-pn) - s1*pn*data.m_mueff;
                
                // reset slip direction
                data.m_s1 = s1;
            }
            else
            {
                t = vec3d(0,0,0);
            }
        }
    }
    // update contact tractions
    else {
        data.m_bstick = false;
        
        if (pme)
        {
            // assume stick and calculate traction
            if (pmep)
            {
                // calculate current global position of the integration point
                vec3d xo = ss.Local2Global(se, n);
                
                // calculate current global position of the previous intersection point
                vec3d xt = ms.Local2Global(*pmep, data.m_rsp[0], data.m_rsp[1]);
                
                // calculate vector gap
                vec3d dg = xt - xo;
                
                // calculate trial stick traction, normal component, shear component
                t = Lt + dg*eps;
                tn = t*nu;
                ts = (t - nu*tn).norm();
                
                // calculate effective friction coefficient
                if (tn != 0)
                {
                    data.m_fls = m_bsmfls ? fls : p/(-tn);
                    mueff = m_mu*(1.0-(1.0-m_phi)*data.m_fls);
                    mueff = MBRACKET(mueff);
                }
                
                // check if stick
                if ( (tn < 0) && (ts < fabs(tn*mueff)) )
                {
                    // set boolean flag for stick
                    data.m_bstick = true;
                    
                    // contact pressure
                    pn = MBRACKET(-tn);
                    
                    // calculate effective friction coefficient
                    if (pn > 0) {
                        data.m_mueff = ts/pn;
                        data.m_fls = m_bsmfls ? fls : p/pn;
                    }

                    // store the previous values as the current
                    data.m_pme = data.m_pmep;
                    data.m_rs = data.m_rsp;
                    
                    // recalculate gaps
                    data.m_dg = dg;
                    
                    // recalculate pressure gap
                    bool mporo = ms.m_poro[pme->m_lid];
                    if (sporo && mporo)
                    {
                        double pm[FEElement::MAX_NODES];
                        for (int k=0; k<pme->Nodes(); ++k) pm[k] = m.Node(pme->m_node[k]).get(m_dofP);
                        double p2 = pme->eval(pm, data.m_rs[0], data.m_rs[1]);
                        data.m_pg = p - p2;
                    }
                    
                }
                else
                {
                    // recalculate contact pressure for slip
                    pn = MBRACKET(Lm + eps*g);
                    
                    if (pn != 0)
                    {
                        
                        double dh = 0;
                        
                        // slip direction
                        s1 = SlipTangent(ss, nel, n, ms, dh, dr);
                        
                        // calculate effective friction coefficient
                        data.m_fls = m_bsmfls ? fls : p/pn;
                        data.m_mueff = m_mu*(1.0-(1.0-m_phi)*data.m_fls);
                        data.m_mueff = MBRACKET(data.m_mueff);

                        // total traction
                        t = nu*(-pn) - s1*pn*data.m_mueff;
                        
                        // reset slip direction
                        data.m_s1 = s1;
                        data.m_bstick = false;
                    }
                    else
                    {
                        t = vec3d(0,0,0);
                    }
                    
                }
            }
            else
            {
                // assume slip upon first contact
                // calculate contact pressure for slip
                pn = MBRACKET(Lm + eps*g);
                
                if (pn != 0)
                {
                    
                    double dh = 0;
                    
                    // slip direction
                    s1 = SlipTangent(ss, nel, n, ms, dh, dr);
                    
                    // calculate effective friction coefficient
                    data.m_fls = m_bsmfls ? fls : p/pn;
                    data.m_mueff = m_mu*(1.0-(1.0-m_phi)*data.m_fls);
                    data.m_mueff = MBRACKET(data.m_mueff);

                    // total traction
                    t = nu*(-pn) - s1*pn*data.m_mueff;
                    
                    // reset slip direction
                    data.m_s1 = s1;
                    data.m_bstick = false;
                }
            }
        }
    }
    
    return t;
    
}

//-----------------------------------------------------------------------------
void FESlidingInterfaceBiphasicMixed::LoadVector(FEGlobalVector& R, const FETimeInfo& tp)
{
	// we will also calculate net contact forces, so zero them here
	m_ss.m_Ft = vec3d(0, 0, 0);
	m_ms.m_Ft = vec3d(0, 0, 0);

	// loop over the nr of passes
	int npass = (m_btwo_pass ? 2 : 1);
	for (int np = 0; np < npass; ++np)
	{
		// get primary and secondary surface
		FESlidingSurfaceBiphasicMixed& ss = (np == 0 ? m_ss : m_ms);
		FESlidingSurfaceBiphasicMixed& ms = (np == 0 ? m_ms : m_ss);

		// assemble the load vector for this pass
		LoadVector(ss, ms, R, tp);
	}
}

//-----------------------------------------------------------------------------
void FESlidingInterfaceBiphasicMixed::LoadVector(FESlidingSurfaceBiphasicMixed& ss, FESlidingSurfaceBiphasicMixed& ms, FEGlobalVector& R, const FETimeInfo& tp)
{
	DOFS& dofs = GetFEModel()->GetDOFS();
	int degree_d = dofs.GetVariableInterpolationOrder(ss.m_varU);
	int degree_p = dofs.GetVariableInterpolationOrder(ss.m_varP);

    const int MN = FEElement::MAX_NODES;
    
    vector<int> sLM, mLM, LM, en;
    vector<double> fe;
    double detJ[MN], w[MN], *Hs, Hm[MN], Hmp[MN];
    double N[4*MN*2];
    
    // need to multiply biphasic force entries by the timestep
    double dt = tp.timeIncrement;
    
    // loop over all primary surface elements
    for (int i=0; i<ss.Elements(); ++i)
    {
        // get the surface element
        FESurfaceElement& se = ss.Element(i);
            
		// flag indicating that primary element is poro
        bool sporo = ss.m_poro[i];
            
        // get the nr of nodes and integration points
        int nseln = se.Nodes();
        int nint = se.GaussPoints();
            
        // copy the LM vector; we'll need it later
        ss.UnpackLM(se, sLM);
            
        // we calculate all the metrics we need before we
        // calculate the nodal forces
        for (int j=0; j<nint; ++j)
        {
            // get the base vectors
            vec3d g[2];
            ss.CoBaseVectors(se, j, g);
                
            // jacobians: J = |g0xg1|
            detJ[j] = (g[0] ^ g[1]).norm();
                
            // integration weights
            w[j] = se.GaussWeights()[j];
        }
            
        // loop over all integration points
        // note that we are integrating over the current surface
        for (int j=0; j<nint; ++j)
        {
            // get the integration point data
			FESlidingSurfaceBiphasicMixed::Data& pt = static_cast<FESlidingSurfaceBiphasicMixed::Data&>(*se.GetMaterialPoint(j));

            // calculate contact pressure and account for stick
            double pn;
            vec3d t = ContactTraction(ss, i, j, ms, pn);
                
            // get the secondary element
            FESurfaceElement* pme = pt.m_pme;
                
            if (pme)
            {
                // get the secondary element
                FESurfaceElement& me = *pme;
                    
				// get the secondary element poro status
                bool mporo = ms.m_poro[pme->m_lid];
                    
                // get the nr of secondary element nodes
                int nmeln = me.Nodes();
                    
                // copy LM vector
                ms.UnpackLM(me, mLM);
                    
                // calculate degrees of freedom
                int ndof = 3*(nseln + nmeln);
                    
                // build the LM vector
                LM.resize(ndof);
                for (int k=0; k<nseln; ++k)
                {
                    LM[3*k  ] = sLM[3*k  ];
                    LM[3*k+1] = sLM[3*k+1];
                    LM[3*k+2] = sLM[3*k+2];
                }
                    
                for (int k=0; k<nmeln; ++k)
                {
                    LM[3*(k+nseln)  ] = mLM[3*k  ];
                    LM[3*(k+nseln)+1] = mLM[3*k+1];
                    LM[3*(k+nseln)+2] = mLM[3*k+2];
                }
                    
                // build the en vector
                en.resize(nseln+nmeln);
                for (int k=0; k<nseln; ++k) en[k      ] = se.m_node[k];
                for (int k=0; k<nmeln; ++k) en[k+nseln] = me.m_node[k];
                    
                // get primary element shape functions
                Hs = se.H(j);
                    
                // get secondary element shape functions
                double r = pt.m_rs[0];
                double s = pt.m_rs[1];
                me.shape_fnc(Hm, r, s);
                    
                if (pn > 0) {
                        
                    // calculate the force vector
                    fe.resize(ndof);
                    zero(fe);
                        
                    for (int k=0; k<nseln; ++k)
                    {
                        N[3*k  ] = Hs[k]*t.x;
                        N[3*k+1] = Hs[k]*t.y;
                        N[3*k+2] = Hs[k]*t.z;
                    }
                        
                    for (int k=0; k<nmeln; ++k)
                    {
                        N[3*(k+nseln)  ] = -Hm[k]*t.x;
                        N[3*(k+nseln)+1] = -Hm[k]*t.y;
                        N[3*(k+nseln)+2] = -Hm[k]*t.z;
                    }
                        
                    for (int k=0; k<ndof; ++k) fe[k] += N[k]*detJ[j]*w[j];
                        
                    // calculate contact forces
                    for (int k=0; k<nseln; ++k)
                        ss.m_Ft += vec3d(fe[3*k], fe[3*k+1], fe[3*k+2]);
                        
                    for (int k = 0; k<nmeln; ++k)
                        ms.m_Ft += vec3d(fe[3*(k+nseln)], fe[3*(k+nseln)+1], fe[3*(k+nseln)+2]);
                        
                    // assemble the global residual
                    R.Assemble(en, LM, fe);
                        
                    // do the biphasic stuff
                    if (sporo && mporo)
                    {
						// get the pressure dofs for each element
						int nspdof = se.ShapeFunctions(degree_p);
						int nmpdof = me.ShapeFunctions(degree_p);
						int npdof = nspdof + nmpdof;

						// evaluate shape functions
						double* Hsp = se.H(degree_p, j);
						me.shape_fnc(degree_p, Hmp, pt.m_rs[0], pt.m_rs[1]);
                            
                        // calculate the flow rate
                        double epsp = m_epsp*pt.m_epsp;
                            
                        double wn = pt.m_Lmp + epsp*pt.m_pg;
                            
                        // fill the LM
                        LM.resize(npdof, -1);
                        for (int k=0; k<nspdof; ++k) LM[k         ] = sLM[3*nseln+k];
                        for (int k=0; k<nmpdof; ++k) LM[k + nspdof] = mLM[3*nmeln+k];
                            
                        // fill the force array
                        fe.resize(npdof);
                        zero(fe);
                        for (int k=0; k<nspdof; ++k) N[k         ] =  Hsp[k];
                        for (int k=0; k<nmpdof; ++k) N[k + nspdof] = -Hmp[k];
                            
                        for (int k=0; k<npdof; ++k) fe[k] += dt*wn*N[k]*detJ[j]*w[j];
                            
						// build the en vector
						en.resize(npdof);
						for (int k = 0; k<nspdof; ++k) en[k         ] = se.m_node[k];
						for (int k = 0; k<nmpdof; ++k) en[k + nspdof] = me.m_node[k];

                        // assemble residual
                        R.Assemble(en, LM, fe);
                    }
                }
            }
        }
    }
}

//-----------------------------------------------------------------------------
void FESlidingInterfaceBiphasicMixed::StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp)
{
	// do single- or two-pass
	int npass = (m_btwo_pass ? 2 : 1);
	for (int np = 0; np < npass; ++np)
	{
		// get the primary and secondary surface
		FESlidingSurfaceBiphasicMixed& ss = (np == 0 ? m_ss : m_ms);
		FESlidingSurfaceBiphasicMixed& ms = (np == 0 ? m_ms : m_ss);

		// assemble stiffness matrix for this pass
		StiffnessMatrix(ss, ms, LS, tp);
	}
}

//-----------------------------------------------------------------------------
//! calculate contact stiffness
void FESlidingInterfaceBiphasicMixed::StiffnessMatrix(FESlidingSurfaceBiphasicMixed& ss, FESlidingSurfaceBiphasicMixed& ms, FELinearSystem& LS, const FETimeInfo& tp)
{
	DOFS& dofs = GetFEModel()->GetDOFS();
	int degree_d = dofs.GetVariableInterpolationOrder(ss.m_varU);
	int degree_p = dofs.GetVariableInterpolationOrder(ss.m_varP);

    // see how many reformations we've had to do so far
    int nref = LS.GetSolver()->m_nref;
    
    const int MN = FEElement::MAX_NODES;
    
    double detJ[MN], w[MN], *Hs, Hm[MN], Hmp[MN];
    double N[4*MN*2], H[4*MN*2];
    vector<int> sLM, mLM, LM, en;
    FEElementMatrix ke;
    
    FEModel& fem = *GetFEModel();
    
    double psf = GetPenaltyScaleFactor();
    
    FEMesh& mesh = *ms.GetMesh();
        
    // loop over all primary surface elements
    for (int i=0; i<ss.Elements(); ++i)
    {
        // get the next element
        FESurfaceElement& se = ss.Element(i);
            
		// primary element's poro status
        bool sporo = ss.m_poro[i];
            
        // get nr of nodes, integration points, pressure dofs
        int nseln = se.Nodes();
        int nint = se.GaussPoints();
		int nspdof = se.ShapeFunctions(degree_p);

        // nodal pressures of primary element
        double pn[MN] = {0};
        if (sporo) {
            for (int j=0; j<nspdof; ++j) pn[j] = ss.GetMesh()->Node(se.m_node[j]).get(m_dofP);
        }
            
        // copy the LM vector
        ss.UnpackLM(se, sLM);

        // we calculate all the metrics we need before we
        // calculate the nodal forces
        for (int j=0; j<nint; ++j)
        {
            // get the base vectors
            vec3d g[2];
            ss.CoBaseVectors(se, j, g);
                
            // jacobians: J = |g0xg1|
            detJ[j] = (g[0] ^ g[1]).norm();
                
            // integration weights
            w[j] = se.GaussWeights()[j];
        }
            
        // loop over all integration points
        for (int j=0; j<nint; ++j)
        {
            // get integration point data
			FESlidingSurfaceBiphasicMixed::Data& pt = static_cast<FESlidingSurfaceBiphasicMixed::Data&>(*se.GetMaterialPoint(j));

            // calculate contact pressure and account for stick
            double pn;
            vec3d t = ContactTraction(ss, i, j, ms, pn);
                
            // get the secondary element
            FESurfaceElement* pme = pt.m_pme;
                
            if (pme)
            {
				// get secondary element
                FESurfaceElement& me = *pme;
                    
				// get secondary element's poro status
                bool mporo = ms.m_poro[pme->m_lid];
                    
                // get the nr of secondary nodes
                int nmeln = me.Nodes();

				// get secondary pressure dofs
				int nmpdof = me.ShapeFunctions(degree_p);
                    
                // nodal pressures of secondary nodes
                double pm[MN] = {0};
                for (int k=0; k<nmpdof; ++k) pm[k] = ms.GetMesh()->Node(me.m_node[k]).get(m_dofP);
                    
                // copy the LM vector
                ms.UnpackLM(me, mLM);
                    
                // calculate degrees of freedom
                int ndpn;	// number of dofs per node
                int ndof;	// number of dofs in stiffness matrix
                    
                if (sporo && mporo) {
                    // calculate degrees of freedom for biphasic-on-biphasic contact
                    ndpn = 4;
                    ndof = ndpn*(nseln+nmeln);
                        
                    // build the LM vector
                    LM.resize(ndof);
                        
                    for (int k=0; k<nseln; ++k)
                    {
                        LM[4*k  ] = sLM[3*k  ];			// x-dof
                        LM[4*k+1] = sLM[3*k+1];			// y-dof
                        LM[4*k+2] = sLM[3*k+2];			// z-dof
                        LM[4*k+3] = sLM[3*nseln+k];		// p-dof
                    }
                    for (int k=0; k<nmeln; ++k)
                    {
                        LM[4*(k+nseln)  ] = mLM[3*k  ];			// x-dof
                        LM[4*(k+nseln)+1] = mLM[3*k+1];			// y-dof
                        LM[4*(k+nseln)+2] = mLM[3*k+2];			// z-dof
                        LM[4*(k+nseln)+3] = mLM[3*nmeln+k];		// p-dof
                    }
                }
                else {
                    // calculate degrees of freedom for biphasic-on-elastic or elastic-on-elastic contact
                    ndpn = 3;
                    ndof = ndpn*(nseln + nmeln);
                        
                    // build the LM vector
                    LM.resize(ndof);
                        
                    for (int k=0; k<nseln; ++k)
                    {
                        LM[3*k  ] = sLM[3*k  ];
                        LM[3*k+1] = sLM[3*k+1];
                        LM[3*k+2] = sLM[3*k+2];
                    }
                        
                    for (int k=0; k<nmeln; ++k)
                    {
                        LM[3*(k+nseln)  ] = mLM[3*k  ];
                        LM[3*(k+nseln)+1] = mLM[3*k+1];
                        LM[3*(k+nseln)+2] = mLM[3*k+2];
                    }
                }
                    
                // build the en vector
                en.resize(nseln+nmeln);
                for (int k=0; k<nseln; ++k) en[k      ] = se.m_node[k];
                for (int k=0; k<nmeln; ++k) en[k+nseln] = me.m_node[k];
                    
                // primary shape functions
                Hs = se.H(j);

				// get primary pressure shape functions
				double* Hsp = se.H(degree_p, j);

                // secondary shape functions
                double r = pt.m_rs[0];
                double s = pt.m_rs[1];
                me.shape_fnc(Hm, r, s);

				// get secondary pressure shape functions
				me.shape_fnc(degree_p, Hmp, r, s);
                    
                // get primary normal vector
                vec3d nu = pt.m_nu;
                    
                // gap function
                double g = pt.m_gap;
                    
                // penalty
                double eps = m_epsn*pt.m_epsn*psf;

                // only evaluate stiffness matrix if contact traction is non-zero
                if (pn > 0)
                {
                    // if stick
                    if (pt.m_bstick)
                    {
                        double dtn = eps;
                            
                        // create the stiffness matrix
                        ke.resize(ndof, ndof); ke.zero();
                            
                        // evaluate basis vectors on primary surface
                        vec3d gscov[2];
                        ss.CoBaseVectors(se, j, gscov);
                            
                        // identity tensor
                        mat3d I = mat3dd(1);
                            
                        // evaluate Mc and Ac and combine them into As
                        double* Gsr = se.Gr(j);
                        double* Gss = se.Gs(j);
                        mat3d Ac[MN], As[MN];
                        mat3d gscovh[2];
                        gscovh[0].skew(gscov[0]); gscovh[1].skew(gscov[1]);
                        for (int k=0; k<nseln; ++k) {
                            Ac[k] = (gscovh[1]*Gsr[k] - gscovh[0]*Gss[k])/detJ[j];
                            As[k] = t & (Ac[k]*nu);
                        }
                            
                        // --- S O L I D - S O L I D   C O N T A C T ---
                            
                        // a. I-term
                        //------------------------------------
                            
                        for (int k=0; k<nseln; ++k) N[k      ] =  Hs[k];
                        for (int k=0; k<nmeln; ++k) N[k+nseln] = -Hm[k];
                            
                        double tmp = dtn*detJ[j]*w[j];
                        for (int l=0; l<nseln+nmeln; ++l)
                        {
                            for (int k=0; k<nseln+nmeln; ++k)
                            {
                                ke[k*ndpn  ][l*ndpn  ] -= -tmp*N[k]*N[l]*I[0][0];
                                ke[k*ndpn  ][l*ndpn+1] -= -tmp*N[k]*N[l]*I[0][1];
                                ke[k*ndpn  ][l*ndpn+2] -= -tmp*N[k]*N[l]*I[0][2];
                                    
                                ke[k*ndpn+1][l*ndpn  ] -= -tmp*N[k]*N[l]*I[1][0];
                                ke[k*ndpn+1][l*ndpn+1] -= -tmp*N[k]*N[l]*I[1][1];
                                ke[k*ndpn+1][l*ndpn+2] -= -tmp*N[k]*N[l]*I[1][2];
                                    
                                ke[k*ndpn+2][l*ndpn  ] -= -tmp*N[k]*N[l]*I[2][0];
                                ke[k*ndpn+2][l*ndpn+1] -= -tmp*N[k]*N[l]*I[2][1];
                                ke[k*ndpn+2][l*ndpn+2] -= -tmp*N[k]*N[l]*I[2][2];
                            }
                        }
                            
                        // b. A-term
                        //-------------------------------------
                            
                        tmp = detJ[j]*w[j];
                        // non-symmetric
                        for (int l=0; l<nseln; ++l)
                        {
                            for (int k=0; k<nseln+nmeln; ++k)
                            {
                                ke[k*ndpn  ][l*ndpn  ] -= tmp*N[k]*As[l][0][0];
                                ke[k*ndpn  ][l*ndpn+1] -= tmp*N[k]*As[l][0][1];
                                ke[k*ndpn  ][l*ndpn+2] -= tmp*N[k]*As[l][0][2];
                                    
                                ke[k*ndpn+1][l*ndpn  ] -= tmp*N[k]*As[l][1][0];
                                ke[k*ndpn+1][l*ndpn+1] -= tmp*N[k]*As[l][1][1];
                                ke[k*ndpn+1][l*ndpn+2] -= tmp*N[k]*As[l][1][2];
                                    
                                ke[k*ndpn+2][l*ndpn  ] -= tmp*N[k]*As[l][2][0];
                                ke[k*ndpn+2][l*ndpn+1] -= tmp*N[k]*As[l][2][1];
                                ke[k*ndpn+2][l*ndpn+2] -= tmp*N[k]*As[l][2][2];
                            }
                        }
                            
                        // --- B I P H A S I C   S T I F F N E S S ---
                        if (sporo && mporo)
                        {
                            // need to multiply biphasic stiffness entries by the timestep
                            double dt = fem.GetTime().timeIncrement;
                                
                            double tmp = dt*w[j]*detJ[j];
                                
                            double epsp = m_epsp*pt.m_epsp*psf;
                                
                            double wn = pt.m_Lmp + epsp*pt.m_pg;
                                
                            // --- S O L I D - P R E S S U R E   C O N T A C T ---
                                
                            // b. A-term
                            //-------------------------------------
                                
                            for (int l=0; l<nseln; ++l) {
                                vec3d Acn = Ac[l]*nu;
                                for (int k=0; k<nseln+nmeln; ++k)
                                {
                                    ke[4*k + 3][4*l  ] -= tmp*wn*N[k]*Acn.x;
                                    ke[4*k + 3][4*l+1] -= tmp*wn*N[k]*Acn.y;
                                    ke[4*k + 3][4*l+2] -= tmp*wn*N[k]*Acn.z;
                                }
                            }
                                
                            // --- P R E S S U R E - P R E S S U R E   C O N T A C T ---
                                
                            // calculate the N-vector
                            for (int k=0; k<nseln; ++k)
                            {
                                N[ndpn*k  ] = 0;
                                N[ndpn*k+1] = 0;
                                N[ndpn*k+2] = 0;
                                N[ndpn*k+3] = Hs[k];
                            }
                                
                            for (int k=0; k<nmeln; ++k)
                            {
                                N[ndpn*(k+nseln)  ] = 0;
                                N[ndpn*(k+nseln)+1] = 0;
                                N[ndpn*(k+nseln)+2] = 0;
                                N[ndpn*(k+nseln)+3] = -Hm[k];
                            }
                                
                            for (int k=0; k<ndof; ++k)
                                for (int l=0; l<ndof; ++l) ke[k][l] -= tmp*epsp*N[k]*N[l];
                                
                        }
                        // assemble the global stiffness
						ke.SetNodes(en);
						ke.SetIndices(LM);
						LS.Assemble(ke);
                    }
                    // if slip
                    else
                    {
                        // create the stiffness matrix
                        ke.resize(ndof, ndof); ke.zero();
                            
                        double tn = -pn;
                            
                        // obtain the slip direction s1 and inverse of spatial increment dh
                        double dh = 0, hd = 0;
                        vec3d dr(0,0,0);
                        vec3d s1 = SlipTangent(ss, i, j, ms, dh, dr);
                            
                        if (dh != 0)
                        {
                            hd = 1.0 / dh;
                        }
                            
                        // evaluate basis vectors on both surfaces
                        vec3d gscov[2], gmcov[2];
                        ss.CoBaseVectors(se, j, gscov);
                        ms.CoBaseVectors(me, r, s, gmcov);
                        mat2d A;
                        A[0][0] = gscov[0]*gmcov[0]; A[0][1] = gscov[0]*gmcov[1];
                        A[1][0] = gscov[1]*gmcov[0]; A[1][1] = gscov[1]*gmcov[1];
                        mat2d a = A.inverse();
                            
                        // evaluate covariant basis vectors on primary surface at previous time step
                        vec3d gscovp[2];
                        ss.CoBaseVectorsP(se, j, gscovp);
                            
                        // calculate delta gscov
                        vec3d dgscov[2];
                        dgscov[0] = gscov[0] - gscovp[0];
                        dgscov[1] = gscov[1] - gscovp[1];
                            
                        // evaluate approximate contravariant basis vectors when gap != 0
                        vec3d gscnt[2], gmcnt[2];
                        gmcnt[0] = gscov[0]*a[0][0] + gscov[1]*a[0][1];
                        gmcnt[1] = gscov[0]*a[1][0] + gscov[1]*a[1][1];
                        gscnt[0] = gmcov[0]*a[0][0] + gmcov[1]*a[1][0];
                        gscnt[1] = gmcov[0]*a[0][1] + gmcov[1]*a[1][1];
                            
                        // evaluate N and S tensors and approximations when gap != 0
                        mat3ds N1 = dyad(nu);
                        mat3d Nh1 = mat3dd(1) - (nu & nu);
                        mat3d Nb1 = mat3dd(1) - (gscov[0] & gscnt[0]) - (gscov[1] & gscnt[1]);
                        mat3d Nt1 = nu & (Nb1*nu);
                        mat3d S1 = s1 & nu;
                        mat3d Sh1 = (mat3dd(1) - (s1 & s1))*hd;
                        mat3d Sb1 = s1 & (Nb1*nu);
                            
                        // evaluate m, c, Mg, and R
                        // evaluate L1 from Mg and R
                        // NOTE: Mg has the 1/detJ included in its definition
                        vec3d m = ((dgscov[0] ^ gscov[1]) + (gscov[0] ^ dgscov[1]));
                        vec3d c = Sh1*Nh1*m*(1/detJ[j]);
                        mat3d Mg = (mat3dd(1)*(nu * m) + (nu & m))*(1/detJ[j]);
                        mat3d B = (c & (Nb1*nu)) - Sh1*Nh1;
                        mat3d R = mat3dd(1)*(nu * dr) + (nu & dr);
                        mat3d L1 = Sh1*((Nh1*Mg - mat3dd(1))*(-g) + R)*Nh1;
                            
                        // evaluate Mc and Ac and combine them into As
                        // evaluate Fc from Ac_bar (Ab)
                        // evaluate Jc as L1*Ac-Fc
                        double* Gsr = se.Gr(j);
                        double* Gss = se.Gs(j);
                        mat3d Ac[MN], As[MN], Pc[MN], Jc[MN];
                        mat3d gscovh[2];
                        mat3d dgscovh[2];
                        gscovh[0].skew(gscov[0]); gscovh[1].skew(gscov[1]);
                        dgscovh[0].skew(dgscov[0]); dgscovh[1].skew(dgscov[1]);
                        for (int k=0; k<nseln; ++k) {
                            vec3d mc = gscnt[0]*Gsr[k] + gscnt[1]*Gss[k];
                            mat3d Mc = nu & mc;
                            Ac[k]    = (gscovh[1]*Gsr[k] - gscovh[0]*Gss[k])/detJ[j];
                            mat3d Ab = (dgscovh[1]*Gsr[k] - dgscovh[0]*Gss[k])/detJ[j];
                            vec3d hcp = (N1*mc + Ac[k]*nu)*pt.m_mueff*(-g);
                            vec3d hcm = (N1*mc*m_mu - Ac[k]*nu*pt.m_mueff);
                            As[k] = (Ac[k] + Mc*N1);
                            mat3d Jc = (L1*Ac[k]) - Sh1*Nh1*Ab*(-g);
                            Pc[k] = (s1 & hcm) + (c & hcp) - Jc*pt.m_mueff;
                        }
                            
                        // evaluate mb and Mb
                        // evaluate s1 dyad mb and combine as Pb
                        double Gmr[MN], Gms[MN];
                        me.shape_deriv(Gmr, Gms, r, s);
                        vec3d mb[MN];
                        mat3d Pb[MN];
                        for (int k=0; k<nmeln; ++k) {
                            mb[k] = gmcnt[0]*Gmr[k] + gmcnt[1]*Gms[k];
                            Pb[k] = ((-nu) & mb[k]) - (s1 & mb[k])*pt.m_mueff;
                        }
                            
                        // evaluate Gbc
                        matrix Gbc(nmeln,nseln);
                        for (int b=0; b<nmeln; ++b) {
                            for (int c=0; c<nseln; ++c) {
                                Gbc(b,c)
                                = (a[0][0]*Gmr[b]*Gsr[c]
                                    + a[0][1]*Gmr[b]*Gss[c]
                                    + a[1][0]*Gms[b]*Gsr[c]
                                    + a[1][1]*Gms[b]*Gss[c])*(-g);
                            }
                        }
                            
                        // define T, Ttb
                        mat3d T = N1 + S1*pt.m_mueff;
                        mat3d Ttb = Nt1 + Sb1*m_mu;

                        // --- S O L I D - S O L I D   C O N T A C T ---
                            
                        // a. NxN-term
                        //------------------------------------
                            
                        for (int k=0; k<nseln; ++k) N[k      ] =  Hs[k];
                        for (int k=0; k<nmeln; ++k) N[k+nseln] = -Hm[k];
                            
                        double tmp = detJ[j]*w[j];
                        for (int l=0; l<nseln+nmeln; ++l)
                        {
                            for (int k=0; k<nseln+nmeln; ++k)
                            {
                                ke[k*ndpn  ][l*ndpn  ] -= -tmp*N[k]*N[l]*(eps*Ttb[0][0] + pt.m_mueff*tn*B[0][0]);
                                ke[k*ndpn  ][l*ndpn+1] -= -tmp*N[k]*N[l]*(eps*Ttb[0][1] + pt.m_mueff*tn*B[0][1]);
                                ke[k*ndpn  ][l*ndpn+2] -= -tmp*N[k]*N[l]*(eps*Ttb[0][2] + pt.m_mueff*tn*B[0][2]);
                                    
                                ke[k*ndpn+1][l*ndpn  ] -= -tmp*N[k]*N[l]*(eps*Ttb[1][0] + pt.m_mueff*tn*B[1][0]);
                                ke[k*ndpn+1][l*ndpn+1] -= -tmp*N[k]*N[l]*(eps*Ttb[1][1] + pt.m_mueff*tn*B[1][1]);
                                ke[k*ndpn+1][l*ndpn+2] -= -tmp*N[k]*N[l]*(eps*Ttb[1][2] + pt.m_mueff*tn*B[1][2]);
                                    
                                ke[k*ndpn+2][l*ndpn  ] -= -tmp*N[k]*N[l]*(eps*Ttb[2][0] + pt.m_mueff*tn*B[2][0]);
                                ke[k*ndpn+2][l*ndpn+1] -= -tmp*N[k]*N[l]*(eps*Ttb[2][1] + pt.m_mueff*tn*B[2][1]);
                                ke[k*ndpn+2][l*ndpn+2] -= -tmp*N[k]*N[l]*(eps*Ttb[2][2] + pt.m_mueff*tn*B[2][2]);
                            }
                        }
                            
                        // b. Na,Nb-term
                        //-------------------------------------
                            
                        tmp = detJ[j]*w[j];
                        // non-symmetric
                        for (int l=0; l<nseln; ++l)
                        {
                            for (int k=0; k<nseln+nmeln; ++k)
                            {
                                ke[k*ndpn  ][l*ndpn  ] -= -tmp*N[k]*(tn*(As[l][0][0] + Pc[l][0][0]));
                                ke[k*ndpn  ][l*ndpn+1] -= -tmp*N[k]*(tn*(As[l][0][1] + Pc[l][0][1]));
                                ke[k*ndpn  ][l*ndpn+2] -= -tmp*N[k]*(tn*(As[l][0][2] + Pc[l][0][2]));
                                    
                                ke[k*ndpn+1][l*ndpn  ] -= -tmp*N[k]*(tn*(As[l][1][0] + Pc[l][1][0]));
                                ke[k*ndpn+1][l*ndpn+1] -= -tmp*N[k]*(tn*(As[l][1][1] + Pc[l][1][1]));
                                ke[k*ndpn+1][l*ndpn+2] -= -tmp*N[k]*(tn*(As[l][1][2] + Pc[l][1][2]));
                                    
                                ke[k*ndpn+2][l*ndpn  ] -= -tmp*N[k]*(tn*(As[l][2][0] + Pc[l][2][0]));
                                ke[k*ndpn+2][l*ndpn+1] -= -tmp*N[k]*(tn*(As[l][2][1] + Pc[l][2][1]));
                                ke[k*ndpn+2][l*ndpn+2] -= -tmp*N[k]*(tn*(As[l][2][2] + Pc[l][2][2]));
                            }
                        }
                            
                        // c. Nc,Nd-term
                        //---------------------------------------
                            
                        tmp = detJ[j]*w[j];
                        // non-symmetric
                        for (int k=0; k<nmeln; ++k)
                        {
                            for (int l=0; l<nseln+nmeln; ++l)
                            {
                                ke[(k+nseln)*ndpn  ][l*ndpn  ] -= tmp*N[l]*tn*Pb[k][0][0];
                                ke[(k+nseln)*ndpn  ][l*ndpn+1] -= tmp*N[l]*tn*Pb[k][0][1];
                                ke[(k+nseln)*ndpn  ][l*ndpn+2] -= tmp*N[l]*tn*Pb[k][0][2];
                                    
                                ke[(k+nseln)*ndpn+1][l*ndpn  ] -= tmp*N[l]*tn*Pb[k][1][0];
                                ke[(k+nseln)*ndpn+1][l*ndpn+1] -= tmp*N[l]*tn*Pb[k][1][1];
                                ke[(k+nseln)*ndpn+1][l*ndpn+2] -= tmp*N[l]*tn*Pb[k][1][2];
                                    
                                ke[(k+nseln)*ndpn+2][l*ndpn  ] -= tmp*N[l]*tn*Pb[k][2][0];
                                ke[(k+nseln)*ndpn+2][l*ndpn+1] -= tmp*N[l]*tn*Pb[k][2][1];
                                ke[(k+nseln)*ndpn+2][l*ndpn+2] -= tmp*N[l]*tn*Pb[k][2][2];
                            }
                        }

                        // c. Gbc-term
                        //---------------------------------------
                            
                        tmp = tn*detJ[j]*w[j];
                        for (int k=0; k<nmeln; ++k)
                        {
                            for (int l=0; l<nseln; ++l)
                            {
                                mat3d gT = T*(Gbc[k][l]*tmp);
                                ke[(k+nseln)*ndpn  ][l*ndpn  ] -= gT[0][0];
                                ke[(k+nseln)*ndpn  ][l*ndpn+1] -= gT[0][1];
                                ke[(k+nseln)*ndpn  ][l*ndpn+2] -= gT[0][2];
                                    
                                ke[(k+nseln)*ndpn+1][l*ndpn  ] -= gT[1][0];
                                ke[(k+nseln)*ndpn+1][l*ndpn+1] -= gT[1][1];
                                ke[(k+nseln)*ndpn+1][l*ndpn+2] -= gT[1][2];
                                    
                                ke[(k+nseln)*ndpn+2][l*ndpn  ] -= gT[2][0];
                                ke[(k+nseln)*ndpn+2][l*ndpn+1] -= gT[2][1];
                                ke[(k+nseln)*ndpn+2][l*ndpn+2] -= gT[2][2];
                            }
                        }

                        // --- B I P H A S I C   S T I F F N E S S ---
                        if (sporo && mporo)
                        {
                            // need to multiply biphasic stiffness entries by the timestep
                            double dt = fem.GetTime().timeIncrement;
                                
                            double dpr = 0, dps = 0;
                            dpr = me.eval_deriv1(degree_p, pm, r, s);
                            dps = me.eval_deriv2(degree_p, pm, r, s);
                                
                            vec3d q2 = gmcnt[0]*dpr + gmcnt[1]*dps;
                                
                            // evaluate gc
                            vector<double> gc(nseln);
                            for (int k=0; k<nseln; ++k) {
                                gc[k]
                                = (a[0][0]*dpr*Gsr[k]
                                    + a[0][1]*dpr*Gss[k]
                                    + a[1][0]*dps*Gsr[k]
                                    + a[1][1]*dps*Gss[k])*(-g);
                            }
                                
                            tmp = dt*w[j]*detJ[j];
                                
                            double epsp = m_epsp*pt.m_epsp*psf;
                                
							// pressure shape functions
							double* Hsp = se.H(degree_p, j);
							me.shape_fnc(degree_p, Hmp, r, s);
							for (int i = 0; i < nseln + nmeln; ++i) H[i] = 0;
							for (int i = 0; i < nspdof; ++i) H[i        ] = Hsp[i];
							for (int i = 0; i < nmpdof; ++i) H[i + nseln] = Hmp[i];
                                
                            // --- S O L I D - P R E S S U R E   C O N T A C T ---
                                
                            // a. q-term
                            //-------------------------------------
                            for (int k=0; k<nseln+nmeln; ++k)
                                for (int l=0; l<nseln+nmeln; ++l)
                                {
                                    ke[4*k + 3][4*l  ] += tmp*epsp*H[k]*N[l]*q2.x;
                                    ke[4*k + 3][4*l+1] += tmp*epsp*H[k]*N[l]*q2.y;
                                    ke[4*k + 3][4*l+2] += tmp*epsp*H[k]*N[l]*q2.z;
                                }
                                
                            double wn = pt.m_Lmp + epsp*pt.m_pg;
                                
                            // b. A-term
                            //-------------------------------------
                                
                            for (int l=0; l<nseln; ++l) {
                                vec3d Acn = Ac[l]*nu;
                                for (int k=0; k<nseln+nmeln; ++k)
                                {
                                    ke[4*k + 3][4*l  ] -= tmp*wn*H[k]*Acn.x;
                                    ke[4*k + 3][4*l+1] -= tmp*wn*H[k]*Acn.y;
                                    ke[4*k + 3][4*l+2] -= tmp*wn*H[k]*Acn.z;
                                }
                            }

                            // c. s-term (Frictional term)
                            //-------------------------------------
                                
							vec3d q = s1*m_mu*(1.0 - m_phi);

                            for (int l=0; l<nseln; ++l) {
                                for (int k=0; k<nseln+nmeln; ++k)
                                {
                                    ke[4*k    ][4*l+3] -= tmp*H[k]*N[l]*q.x;
                                    ke[4*k + 1][4*l+3] -= tmp*H[k]*N[l]*q.y;
                                    ke[4*k + 2][4*l+3] -= tmp*H[k]*N[l]*q.z;
                                }
                            }

                            // d. m-term
                            //---------------------------------------
                                
							vec3d mbp[MN];
							double Gpr[MN], Gps[MN];
							me.shape_deriv(degree_p, Gpr, Gps, r, s);
							for (int k = 0; k < nmeln; ++k) mbp[k] = vec3d(0, 0, 0);
							for (int k = 0; k<nmpdof; ++k)
							{
								mbp[k] = gmcnt[0] * Gpr[k] + gmcnt[1] * Gps[k];
							}

                            for (int k=0; k<nmpdof; ++k) {
                                for (int l=0; l<nseln+nmeln; ++l)
                                {
                                    ke[4*(k+nseln) + 3][4*l  ] += tmp*wn*N[l]*mbp[k].x;
                                    ke[4*(k+nseln) + 3][4*l+1] += tmp*wn*N[l]*mbp[k].y;
                                    ke[4*(k+nseln) + 3][4*l+2] += tmp*wn*N[l]*mbp[k].z;
                                }
                            }

                            // e. gc-term
                            //-------------------------------------
                            for (int k=0; k<nseln+nmeln; ++k)
                                for (int l=0; l<nseln; ++l)
                                {
                                    ke[4*k + 3][4*l  ] -= tmp*epsp*H[k]*gc[l]*nu.x;
                                    ke[4*k + 3][4*l+1] -= tmp*epsp*H[k]*gc[l]*nu.y;
                                    ke[4*k + 3][4*l+2] -= tmp*epsp*H[k]*gc[l]*nu.z;
                                }
                                
                            // f. Gbc-term (CONVERT!)
                            //---------------------------------------
/*                                
                            for (int k=0; k<nmeln; ++k) {
                                for (int l=0; l<nseln; ++l)
                                {
                                    ke[4*(k+nseln) + 3][4*l  ] -= tmp*wn*Gbc[k][l]*nu.x;
                                    ke[4*(k+nseln) + 3][4*l+1] -= tmp*wn*Gbc[k][l]*nu.y;
                                    ke[4*(k+nseln) + 3][4*l+2] -= tmp*wn*Gbc[k][l]*nu.z;
                                }
                            }
*/
                            // --- P R E S S U R E - P R E S S U R E   C O N T A C T ---
                                
                            // calculate the N-vector
							for (int k = 0; k < ndof; ++k) N[k] = 0.0;
                            for (int k=0; k<nspdof; ++k)
                            {
                                N[ndpn*k+3] = Hsp[k];
                            }
                                
                            for (int k=0; k<nmpdof; ++k)
                            {
                                N[ndpn*(k+nseln)+3] = -Hmp[k];
                            }
                                
                            for (int k=0; k<ndof; ++k)
                                for (int l=0; l<ndof; ++l) ke[k][l] -= tmp*epsp*N[k]*N[l];
                        }
                            
                        // assemble the global stiffness
						ke.SetNodes(en);
						ke.SetIndices(LM);
						LS.Assemble(ke);
                    }
                }
            }
        }
    }
}

//-----------------------------------------------------------------------------
void FESlidingInterfaceBiphasicMixed::UpdateContactPressures()
{
	DOFS& dofs = GetFEModel()->GetDOFS();
	int degree_d = dofs.GetVariableInterpolationOrder(m_ss.m_varU);
	int degree_p = dofs.GetVariableInterpolationOrder(m_ss.m_varP);

    int npass = (m_btwo_pass?2:1);
    const int MN = FEElement::MAX_NODES;
    const int MI = FEElement::MAX_INTPOINTS;
    
    double psf = GetPenaltyScaleFactor();
    
    for (int np=0; np<npass; ++np)
    {
		FESlidingSurfaceBiphasicMixed& ss = (np == 0? m_ss : m_ms);
		FESlidingSurfaceBiphasicMixed& ms = (np == 0? m_ms : m_ss);
        
        // loop over all elements of the primary surface
        for (int n=0; n<ss.Elements(); ++n)
        {
            FESurfaceElement& el = ss.Element(n);
            int nint = el.GaussPoints();
            
            // get the normal tractions at the integration points
            for (int i=0; i<nint; ++i)
            {
                // get integration point data
				FESlidingSurfaceBiphasicMixed::Data& sd = static_cast<FESlidingSurfaceBiphasicMixed::Data&>(*el.GetMaterialPoint(i));
				// evaluate traction on primary surface
                double eps = m_epsn*sd.m_epsn*psf;
                if (sd.m_bstick) {
                    // if stick, evaluate total traction
                    sd.m_tr = sd.m_Lmt + sd.m_dg*eps;
                    // then derive normal component
                    sd.m_Ln = -sd.m_tr*sd.m_nu;
                }
                else {
                    // if slip, evaluate normal traction
                    double Ln = sd.m_Lmd + eps*sd.m_gap;
                    sd.m_Ln = MBRACKET(Ln);
                    // then derive total traction
                    sd.m_tr = -(sd.m_nu*sd.m_Ln + sd.m_s1*sd.m_Ln*sd.m_mueff);
                    
                }
                
                FESurfaceElement* pme = sd.m_pme;
                
                if (m_btwo_pass && pme)
                {
                    // get secondary element data
                    int mint = pme->GaussPoints();
                    double pi[MI];
                    vec3d ti[MI];
                    for (int j=0; j<mint; ++j)
                    {
						FESlidingSurfaceBiphasicMixed::Data& md = static_cast<FESlidingSurfaceBiphasicMixed::Data&>(*pme->GetMaterialPoint(j));

                        // evaluate traction on secondary surface
                        double eps = m_epsn*md.m_epsn*psf;
                        if (md.m_bstick) {
                            // if stick, evaluate total traction
                            ti[j] = md.m_Lmt + md.m_dg*eps;
                            // then derive normal component
                            pi[j] = -ti[j]*md.m_nu;
                        }
                        else {
                            // if slip, evaluate normal traction
                            double Ln = md.m_Lmd + eps*md.m_gap;
                            pi[j] = MBRACKET(Ln);
                            // then derive total traction
                            ti[j] = -(md.m_nu*pi[j] + md.m_s1*md.m_mueff*pi[j]);
                        }
                    }
                    // project the data to the nodes
                    double pn[MN];
                    vec3d tn[MN];
                    pme->FEElement::project_to_nodes(pi, pn);
                    pme->project_to_nodes(ti, tn);
                    // now evaluate the traction at the intersection point
                    double Ln = pme->eval(degree_p, pn, sd.m_rs[0], sd.m_rs[1]);
                    vec3d trac = pme->eval(tn, sd.m_rs[0], sd.m_rs[1]);
                    sd.m_Ln += MBRACKET(Ln);
                    // tractions on primary-secondary are opposite, so subtract
                    sd.m_tr -= trac;
                }
            }
        }
        ss.EvaluateNodalContactPressures();
        ss.EvaluateNodalContactTractions();
    }
}

//-----------------------------------------------------------------------------
bool FESlidingInterfaceBiphasicMixed::Augment(int naug, const FETimeInfo& tp)
{
    // make sure we need to augment
	if (m_laugon != 1) return true;

    int i;
    double Ln, Lp;
    bool bconv = true;
    
    double psf = GetPenaltyScaleFactor();
    
    bool bporo = (m_ss.m_bporo && m_ms.m_bporo);
    int NS = m_ss.Elements();
    int NM = m_ms.Elements();
    
    // --- c a l c u l a t e   i n i t i a l   n o r m s ---
    // a. normal component
    double normL0 = 0, normP = 0, normDP = 0;
    for (int i=0; i<NS; ++i)
    {
		FESurfaceElement& el = m_ss.Element(i);
        for (int j=0; j<el.GaussPoints(); ++j)
        {
			FESlidingSurfaceBiphasicMixed::Data& ds = static_cast<FESlidingSurfaceBiphasicMixed::Data&>(*el.GetMaterialPoint(j));
			if (ds.m_bstick)
                normL0 += ds.m_Lmt*ds.m_Lmt;
            else
                normL0 += ds.m_Lmd*ds.m_Lmd;
        }
    }
    for (int i=0; i<NM; ++i)
    {
		FESurfaceElement& el = m_ms.Element(i);
        for (int j=0; j<el.GaussPoints(); ++j)
        {
			FESlidingSurfaceBiphasicMixed::Data& dm = static_cast<FESlidingSurfaceBiphasicMixed::Data&>(*el.GetMaterialPoint(j));
			if (dm.m_bstick)
                normL0 += dm.m_Lmt*dm.m_Lmt;
            else
                normL0 += dm.m_Lmd*dm.m_Lmd;
        }
    }
    
    // b. gap component
    // (is calculated during update)
    double maxgap = 0;
    double maxpg = 0;
    
    // update Lagrange multipliers
    double normL1 = 0, epsp;
    for (i=0; i<NS; ++i)
    {
		FESurfaceElement& el = m_ss.Element(i);
        vec3d tn[FEElement::MAX_INTPOINTS];
        if (m_bsmaug) m_ss.GetGPSurfaceTraction(i, tn);
        for (int j=0; j<el.GaussPoints(); ++j)
        {
			FESlidingSurfaceBiphasicMixed::Data& ds = static_cast<FESlidingSurfaceBiphasicMixed::Data&>(*el.GetMaterialPoint(j));

            // update Lagrange multipliers on primary surface
            double eps = m_epsn*ds.m_epsn*psf;
            if (ds.m_bstick) {
                // if stick, augment total traction
                if (m_bsmaug) {
                    ds.m_Lmt = tn[j];
                    if (m_btwo_pass) ds.m_Lmt /= 2;
                }
                else {
                    ds.m_Lmt += ds.m_dg*eps;
                }
                // then derive normal component
                ds.m_Lmd = -ds.m_Lmt*ds.m_nu;
                Ln = ds.m_Lmd;
                normL1 += ds.m_Lmt*ds.m_Lmt;
                
                if (Ln > 0) maxgap = max(maxgap, fabs(ds.m_dg.norm()));
            }
            else {
                // if slip, augment normal traction
                if (m_bsmaug) {
                    Ln = -(tn[j]*ds.m_nu);
                    ds.m_Lmd = MBRACKET(Ln);
                    if (m_btwo_pass) ds.m_Lmd /= 2;
                }
                else {
                    Ln = ds.m_Lmd + eps*ds.m_gap;
                    ds.m_Lmd = MBRACKET(Ln);
                }
                // then derive total traction
                double mueff = m_mu*(1.0-(1.0-m_phi)*ds.m_p1/ds.m_Lmd);
                if ( ds.m_Lmd < (1-m_phi)*ds.m_p1 )
                {
                    mueff = 0.0;
                }
                ds.m_Lmt = -(ds.m_nu*ds.m_Lmd + ds.m_s1*ds.m_Lmd*mueff);
                normL1 += ds.m_Lmd*ds.m_Lmd;
                
                if (Ln > 0) maxgap = max(maxgap, fabs(ds.m_gap));
            }
            if (m_ss.m_bporo) {
                Lp = 0;
                if (Ln > 0) {
                    epsp = m_epsp*ds.m_epsp*psf;
                    Lp = ds.m_Lmp + epsp*ds.m_pg;
                    maxpg = max(maxpg,fabs(ds.m_pg));
                    normDP += ds.m_pg*ds.m_pg;
                }
                ds.m_Lmp = Lp;
            }
        }
    }
    
    for (i=0; i<NM; ++i)
    {
		FESurfaceElement& el = m_ms.Element(i);
        vec3d tn[FEElement::MAX_INTPOINTS];
        if (m_bsmaug) m_ms.GetGPSurfaceTraction(i, tn);
        for (int j=0; j<el.GaussPoints(); ++j)
        {
			FESlidingSurfaceBiphasicMixed::Data& dm = static_cast<FESlidingSurfaceBiphasicMixed::Data&>(*el.GetMaterialPoint(j));

            // update Lagrange multipliers on secondary surface
            double eps = m_epsn*dm.m_epsn*psf;
            if (dm.m_bstick) {
                // if stick, augment total traction
                if (m_bsmaug) {
                    dm.m_Lmt = tn[j];
                    if (m_btwo_pass) dm.m_Lmt /= 2;
                }
                else {
                    dm.m_Lmt += dm.m_dg*eps;
                }
                // then derive normal component
                dm.m_Lmd = -dm.m_Lmt*dm.m_nu;
                Ln = dm.m_Lmd;
                normL1 += dm.m_Lmt*dm.m_Lmt;
                
                if (Ln > 0) maxgap = max(maxgap, fabs(dm.m_dg.norm()));
            }
            else {
                // if slip, augment normal traction
                if (m_bsmaug) {
                    Ln = -(tn[j]*dm.m_nu);
                    dm.m_Lmd = MBRACKET(Ln);
                    if (m_btwo_pass) dm.m_Lmd /= 2;
                }
                else {
                    Ln = dm.m_Lmd + eps*dm.m_gap;
                    dm.m_Lmd = MBRACKET(Ln);
                }
                // then derive total traction
                double mueff = m_mu*(1.0-(1.0-m_phi)*dm.m_p1/dm.m_Lmd);
                if ( dm.m_Lmd < (1-m_phi)*dm.m_p1 )
                {
                    mueff = 0.0;
                }
                dm.m_Lmt = -(dm.m_nu*dm.m_Lmd + dm.m_s1*dm.m_Lmd*mueff);
                normL1 += dm.m_Lmd*dm.m_Lmd;
                
                if (Ln > 0) maxgap = max(maxgap, fabs(dm.m_gap));
            }
            
            if (m_ms.m_bporo) {
                Lp = 0;
                if (Ln > 0) {
                    epsp = m_epsp*dm.m_epsp*psf;
                    Lp = dm.m_Lmp + epsp*dm.m_pg;
                    maxpg = max(maxpg,fabs(dm.m_pg));
                    normDP += dm.m_pg*dm.m_pg;
                }
                dm.m_Lmp = Lp;
            }
        }
    }
    
    // normP should be a measure of the fluid pressure at the
    // contact interface.  However, since it could be zero,
    // use an average measure of the contact traction instead.
    normP = normL1;
    
    // calculate relative norms
    double lnorm = (normL1 != 0 ? fabs((normL1 - normL0) / normL1) : fabs(normL1 - normL0));
    double pnorm = (normP != 0 ? (normDP/normP) : normDP);
    
    // check convergence
    if ((m_gtol > 0) && (maxgap > m_gtol)) bconv = false;
    if ((m_ptol > 0) && (bporo && maxpg > m_ptol)) bconv = false;
    
    if ((m_atol > 0) && (lnorm > m_atol)) bconv = false;
    if ((m_atol > 0) && (pnorm > m_atol)) bconv = false;
    
    if (naug < m_naugmin ) bconv = false;
    if (naug >= m_naugmax) bconv = true;
    
    feLog(" sliding interface # %d\n", GetID());
    feLog("                        CURRENT        REQUIRED\n");
    feLog("    D multiplier : %15le", lnorm); if (m_atol > 0) feLog("%15le\n", m_atol); else feLog("       ***\n");
    if (bporo) { feLog("    P gap        : %15le", pnorm); if (m_atol > 0) feLog("%15le\n", m_atol); else feLog("       ***\n"); }
    
    feLog("    maximum gap  : %15le", maxgap);
    if (m_gtol > 0) feLog("%15le\n", m_gtol); else feLog("       ***\n");
    if (bporo) {
        feLog("    maximum pgap : %15le", maxpg);
        if (m_ptol > 0) feLog("%15le\n", m_ptol); else feLog("       ***\n");
    }
    
    ProjectSurface(m_ss, m_ms, true);
    if (m_btwo_pass) ProjectSurface(m_ms, m_ss, true);
    
    m_bfreeze = true;
    
    return bconv;
}

//-----------------------------------------------------------------------------
void FESlidingInterfaceBiphasicMixed::Serialize(DumpStream &ar)
{
    // serialize contact data
    FEContactInterface::Serialize(ar);
    
    // serialize contact surface data
    m_ms.Serialize(ar);
    m_ss.Serialize(ar);

	// serialize element pointers
	SerializeElementPointers(m_ss, m_ms, ar);
	SerializeElementPointers(m_ms, m_ss, ar);
}

//-----------------------------------------------------------------------------

void FESlidingInterfaceBiphasicMixed::MarkFreeDraining()
{
    // Mark all nodes as free-draining.  This needs to be done for ALL
    // contact interfaces prior to executing Update(), where nodes that are
    // in contact are subsequently marked as non free-draining.  This ensures
    // that for surfaces involved in more than one contact interface, nodes
    // that have been marked as non free-draining are not reset to
    // free-draining.
    for (int np=0; np<2; ++np)
    {
		FESlidingSurfaceBiphasicMixed& s = (np == 0? m_ss : m_ms);
        
        if (s.m_bporo) {
            // first, mark all nodes as free-draining (= neg. ID)
            // this is done by setting the dof's equation number
            // to a negative number
            for (int i=0; i<s.Nodes(); ++i)
            {
				FENode& node = s.Node(i);
				int id = node.m_ID[m_dofP];
                if (id >= 0)
                {
                    // mark node as free-draining
                    node.m_ID[m_dofP] = -id-2;
                }
            }
        }
    }
}

//-----------------------------------------------------------------------------
void FESlidingInterfaceBiphasicMixed::SetFreeDraining()
{
    // Set the pressure to zero for the free-draining nodes
    for (int np=0; np<2; ++np)
    {
        FESlidingSurfaceBiphasicMixed& s = (np == 0? m_ss : m_ms);
        
        if (s.m_bporo) {
            // loop over all nodes
            for (int i=0; i<s.Nodes(); ++i)
            {
				FENode& node = s.Node(i);
				if (node.m_ID[m_dofP] < -1)
                {
                    // set the fluid pressure to zero
                    node.set(m_dofP, 0);
                }
            }
        }
    }
}
