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
#include "FEBioFSI.h"
#include <FECore/FEModel.h>
#include <FECore/FEAnalysis.h>
#include <FECore/FENormalProjection.h>
#include <FECore/FELinearSystem.h>
#include <FECore/log.h>

//-----------------------------------------------------------------------------
// Define sliding interface parameters
BEGIN_FECORE_CLASS(FETiedFluidFSI, FEContactInterface)
	ADD_PARAMETER(m_laugon   , "laugon"             )->setLongName("Enforcement method")->setEnums("PENALTY\0AUGLAG\0");
	ADD_PARAMETER(m_atol     , "tolerance"          );
	ADD_PARAMETER(m_gtol     , "gaptol"             );
	ADD_PARAMETER(m_epsn     , "penalty"            );
	ADD_PARAMETER(m_bautopen , "auto_penalty"       );
    ADD_PARAMETER(m_bupdtpen , "update_penalty"     );
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
    m_dg = vec3d(0,0,0);
    m_nu = vec3d(0,0,0);
    m_rs = vec2d(0,0);
    m_Lmd = vec3d(0,0,0);
    m_tr = vec3d(0,0,0);
    m_epsn = 1.0;
}

//-----------------------------------------------------------------------------
void FETiedFluidFSISurface::Data::Serialize(DumpStream& ar)
{
	FEContactMaterialPoint::Serialize(ar);
	ar & m_Gap;
	ar & m_dg;
	ar & m_nu;
	ar & m_rs;
	ar & m_Lmd;
	ar & m_epsn;
	ar & m_tr;
}

void FETiedFluidFSISurface::Data::Init()
{
    FEContactMaterialPoint::Init();
    m_Gap = vec3d(0, 0, 0);
    m_dg = vec3d(0, 0, 0);
    m_nu = vec3d(0, 0, 0);
    m_rs = vec2d(0, 0);
    m_Lmd = vec3d(0, 0, 0);
    m_tr = vec3d(0, 0, 0);
    m_epsn = 1.0;
}

//-----------------------------------------------------------------------------
// FETiedFluidFSISurface
//-----------------------------------------------------------------------------

FETiedFluidFSISurface::FETiedFluidFSISurface(FEModel* pfem) : FEContactSurface(pfem)
{

}

//-----------------------------------------------------------------------------
bool FETiedFluidFSISurface::Init()
{
    // initialize surface data first
    if (FEContactSurface::Init() == false) return false;
    
    // allocate node normals
    m_nn.assign(Nodes(), vec3d(0,0,0));
    
    // initialize nodal force vector
    m_Fn.assign(Nodes(), vec3d(0, 0, 0));

    return true;
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
void FETiedFluidFSISurface::GetContactTraction(int nface, vec3d& pt)
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
vec3d FETiedFluidFSISurface::GetContactForce()
{
    // initialize contact force
    vec3d f(0, 0, 0);

    // loop over all elements of the primary surface
    for (int i = 0; i < m_Fn.size(); ++i) f += m_Fn[i];

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
            double T = data.m_tr.norm2();
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
    m_epsn = 1;
    m_stol = 0.01;
    m_bsymm = true;
    m_srad = 1.0;
    m_gtol = -1;    // we use augmentation tolerance by default
    m_bautopen = false;
    m_bupdtpen = false;
    
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
    
    // get the degrees of freedom
    // TODO: Can this be done in Init, since  there is no error checking
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
}

//-----------------------------------------------------------------------------
bool FETiedFluidFSI::Init()
{
    // initialize surface data
    if (m_ss.Init() == false) return false;
    if (m_ms.Init() == false) return false;
    m_ss.SetShellBottom(m_bshellb);

	// Flip secondary and primary surfaces, if requested.
	// Note that we turn off those flags because otherwise we keep flipping, each time we get here (e.g. in optimization)
	// TODO: Of course, we shouldn't get here more than once. I think we also get through the FEModel::Reset, so I'll have
	//       look into that. 
	if (m_bflips) { m_ss.Invert(); m_bflips = false; }
	if (m_bflipm) { m_ms.Invert(); m_bflipm = false; }
    
    return true;
}

//-----------------------------------------------------------------------------
//! build the matrix profile for use in the stiffness matrix
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

    const int ndpn1 = 3;
    const int ndpn2 = 7;
    vector<int> lm((ndpn1+ndpn2)*FEElement::MAX_NODES);
    
    for (int np=0; np<1; ++np)
    {
        FETiedFluidFSISurface& ss = (np == 0? m_ss : m_ms);
        
        int ni = 0, k, l;
        for (int j=0; j<ss.Elements(); ++j)
        {
            FESurfaceElement& se = ss.Element(j);
            int nint = se.GaussPoints();
            int* sn = &se.m_node[0];
            for (k=0; k<nint; ++k, ++ni)
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
                    
                    for (l=0; l<nseln; ++l)
                    {
                        vector<int>& id = mesh.Node(sn[l]).m_ID;
                        lm[ndpn1*l  ] = id[dof_X];
                        lm[ndpn1*l+1] = id[dof_Y];
                        lm[ndpn1*l+2] = id[dof_Z];
                    }
                    
                    for (l=0; l<nmeln; ++l)
                    {
                        vector<int>& id = mesh.Node(mn[l]).m_ID;
                        lm[ndpn1*(l+nseln)  ] = id[dof_X];
                        lm[ndpn1*(l+nseln)+1] = id[dof_Y];
                        lm[ndpn1*(l+nseln)+2] = id[dof_Z];
                        lm[ndpn1*(l+nseln)+3] = id[dof_WX];
                        lm[ndpn1*(l+nseln)+4] = id[dof_WY];
                        lm[ndpn1*(l+nseln)+5] = id[dof_WZ];
                        lm[ndpn1*(l+nseln)+6] = id[dof_EF];
                    }
                    
                    K.build_add(lm);
                }
            }
        }
    }
}

//-----------------------------------------------------------------------------
void FETiedFluidFSI::UpdateAutoPenalty()
{
    // calculate the penalty
    if (m_bautopen)
    {
        CalcAutoPenalty(m_ss);
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
    InitialProjection(m_ss, m_ms);
    
    // get the list of fluid-FSI elements connected to this interface
    FEModel* fem = GetFEModel();
    int NF = m_ms.Elements();
    m_elem.resize(NF);
    m_s.resize(NF, 1);
    for (int j = 0; j < NF; ++j)
    {
        bool bself = false;
        FESurfaceElement& el = m_ms.Element(j);
        // extract the first of two elements on this interface
        m_elem[j] = el.m_elem[0].pe;
        // get its material and check if FluidFSI
        FEMaterial* pm = fem->GetMaterial(m_elem[j]->GetMatID());
        FEFluidFSI* pfsi = dynamic_cast<FEFluidFSI*>(pm);
        if (pfsi) {
            double s = m_ms.FacePointing(el, *m_elem[j]);
            m_s[j] = bself ? -s : s;
            assert(m_s[j]);
        }
        else if (!bself) {
            // extract the second of two elements on this interface
            m_elem[j] = el.m_elem[1].pe;
            pm = fem->GetMaterial(m_elem[j]->GetMatID());
            pfsi = dynamic_cast<FEFluidFSI*>(pm);
            assert(pfsi);
            m_s[j] = m_ms.FacePointing(el, *m_elem[j]);
            assert(m_s[j]);
        }
        else
            assert(false);
    }
}

//-----------------------------------------------------------------------------
void FETiedFluidFSI::CalcAutoPenalty(FETiedFluidFSISurface& s)
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
            FETiedFluidFSISurface::Data& pt = static_cast<FETiedFluidFSISurface::Data&>(*el.GetMaterialPoint(j));
			pt.m_epsn = eps;
        }
    }
}

//-----------------------------------------------------------------------------
// Perform initial projection between tied surfaces in reference configuration
void FETiedFluidFSI::InitialProjection(FETiedFluidFSISurface& ss, FETiedFluidFSISurface& ms)
{
    FESurfaceElement* pme;
    vec3d r, nu;
    double rs[2];
    
    // initialize projection data
    FENormalProjection np(ms);
    np.SetTolerance(m_stol);
    np.SetSearchRadius(m_srad);
    np.Init();

	// let's count the number of contact pairs we find.
	int contacts = 0;
    
    // loop over all integration points
    int n = 0;
    for (int i=0; i<ss.Elements(); ++i)
    {
        FESurfaceElement& el = ss.Element(i);
        
        int nint = el.GaussPoints();
        
        for (int j=0; j<nint; ++j, ++n)
        {
            // calculate the global position of the integration point
            r = ss.Local2Global(el, j);
            
            // calculate the normal at this integration point
            nu = ss.SurfaceNormal(el, j);
            
            // find the intersection point with the secondary surface
            pme = np.Project2(r, nu, rs);
            
            FETiedFluidFSISurface::Data& pt = static_cast<FETiedFluidFSISurface::Data&>(*el.GetMaterialPoint(j));
			pt.m_pme = pme;
            pt.m_rs[0] = rs[0];
            pt.m_rs[1] = rs[1];
            if (pme)
            {
                // the node could potentially be in contact
                // find the global location of the intersection point
                vec3d q = ms.Local2Global(*pme, rs[0], rs[1]);
                
                // calculate the gap function
                pt.m_Gap = q - r;

				contacts++;
            }
            else
            {
                // the node is not in contact
                pt.m_Gap = vec3d(0,0,0);
            }
        }
    }

	// if we found no contact pairs, let's report this since this is probably not the user's intention
	if (contacts == 0)
	{
		std::string name = GetName();
		feLogWarning("No contact pairs found for tied interface \"%s\".\nThis contact interface may not have any effect.", name.c_str());
	}
}

//-----------------------------------------------------------------------------
// Evaluate gap functions for position and fluid pressure
void FETiedFluidFSI::ProjectSurface(FETiedFluidFSISurface& ss, FETiedFluidFSISurface& ms)
{
    FESurfaceElement* pme;
    vec3d r;
    
    // loop over all integration points
    for (int i=0; i<ss.Elements(); ++i)
    {
        FESurfaceElement& el = ss.Element(i);
        
        int nint = el.GaussPoints();
        
        for (int j=0; j<nint; ++j)
        {
            FETiedFluidFSISurface::Data& pt = static_cast<FETiedFluidFSISurface::Data&>(*el.GetMaterialPoint(j));

            // calculate the global position of the integration point
            r = ss.Local2Global(el, j);
            
            // calculate the normal at this integration point
            pt.m_nu = ss.SurfaceNormal(el, j);
            
            // if this node is tied, evaluate gap functions
            pme = pt.m_pme;
            if (pme)
            {
                // find the global location of the intersection point
                vec3d q = ms.Local2Global(*pme, pt.m_rs[0], pt.m_rs[1]);
                
                // calculate the gap function
                vec3d g = q - r;
                pt.m_dg = g - pt.m_Gap;
				pt.m_gap = pt.m_dg.norm();                
            }
            else
            {
                // the node is not tied
                pt.m_dg = vec3d(0,0,0);
				pt.m_gap = 0.0;
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
    
}

//-----------------------------------------------------------------------------
double FETiedFluidFSI::GetFluidDilatation(FESurfaceElement& el, vec2d rs, double alpha)
{
    double ef = 0;
    double* H;
    el.shape_fnc(H, rs.x(), rs.y());
    int neln = el.Nodes();
    for (int j = 0; j < neln; ++j) {
        FENode& node = m_ms.Node(el.m_lnode[j]);
        double ej = node.get(m_dofEF)*alpha + node.get_prev(m_dofEF)*(1.0 - alpha);
        ef += ej*H[j];
    }
    return ef;
}

//-----------------------------------------------------------------------------
mat3ds FETiedFluidFSI::GetFluidStress(FESurfaceElement& el, vec2d rs)
{
    FEElement* e = el.m_elem[0].pe;
    // get the fluid-FSI material of this element
    FEFluidFSI* pfsi = dynamic_cast<FEFluidFSI*>(GetFEModel()->GetMaterial(e->GetMatID()));

    if (e && pfsi) {
        mat3ds si[FEElement::MAX_INTPOINTS];
        mat3ds so[FEElement::MAX_NODES];
        for (int i=0; i<el.GaussPoints(); ++i) {
            FEMaterialPoint* pt = el.GetMaterialPoint(i);
            FEFluidMaterialPoint* ep = pt->ExtractData<FEFluidMaterialPoint>();
            if (ep)
                si[i] = ep->m_sf;
            else
                si[i].zero();
        }
        // project stresses from integration points to nodes
        el.project_to_nodes(si, so);
        // only keep the stresses at the nodes of the contact face
        mat3ds sn[FEElement::MAX_NODES];
        for (int i=0; i<el.Nodes(); ++i)
            sn[i] = so[el.FindNode(el.m_node[i])];
        // evaluate fluid stress at given parametric coordinates
        mat3ds sf; sf.zero();
        double H[FEElement::MAX_INTPOINTS];
        el.shape_fnc(H, rs.x(), rs.y());
        for (int j=0; j<el.Nodes(); ++j) {
            sf += sn[j]*H[j];
        }
        return sf;
    }
    else
        return mat3ds(0);
}

//-----------------------------------------------------------------------------
mat3ds FETiedFluidFSI::GetViscousFluidStress(FESurfaceElement& el, vec2d rs)
{
    FEElement* e = el.m_elem[0].pe;
    // get the fluid-FSI material of this element
    FEFluidFSI* pfsi = dynamic_cast<FEFluidFSI*>(GetFEModel()->GetMaterial(e->GetMatID()));
    FESolidElement* se = dynamic_cast<FESolidElement*>(e);
    
    if (se && pfsi) {
        mat3ds si[FEElement::MAX_INTPOINTS];
        mat3ds so[FEElement::MAX_NODES];
        for (int i=0; i<se->GaussPoints(); ++i) {
            FEMaterialPoint& pt = *se->GetMaterialPoint(i);
            si[i] = pfsi->Fluid()->GetViscous()->Stress(pt);
        }
        // project stresses from integration points to nodes
        se->project_to_nodes(si, so);
        // only keep the stresses at the nodes of the contact face
        mat3ds sn[FEElement::MAX_NODES];
        for (int i=0; i<el.Nodes(); ++i)
            sn[i] = so[se->FindNode(el.m_node[i])];
        // evaluate fluid stress at given parametric coordinates
        mat3ds sv; sv.zero();
        double H[FEElement::MAX_INTPOINTS];
        el.shape_fnc(H, rs.x(), rs.y());
        for (int j=0; j<el.Nodes(); ++j) {
            sv += sn[j]*H[j];
        }
        return sv;
    }
    else
        return mat3ds(0);
}

//-----------------------------------------------------------------------------
void FETiedFluidFSI::LoadVector(FEGlobalVector& R, const FETimeInfo& tp)
{
    vector<int> LM1, LM2, LM, en;
    vector<double> fe;
    const int MI = FEElement::MAX_INTPOINTS;
    double detJ[MI], w[MI], *N1, N2[MI];
    vec3d n1[MI];
    const int MN = FEElement::MAX_NODES;
    vec3d f1[MN], f2[MN];

    // zero nodal forces
    m_ss.m_Fn.assign(m_ss.Nodes(), vec3d(0, 0, 0));
    m_ms.m_Fn.assign(m_ms.Nodes(), vec3d(0, 0, 0));

    // Update auto-penalty if requested
    if (m_bupdtpen && (GetFEModel()->GetCurrentStep()->GetFESolver()->m_niter == 0))
        UpdateAutoPenalty();
    
    // loop over the nr of passes
    int npass = 1;
    for (int np=0; np<npass; ++np)
    {
        // get primary and secondary surface
        FETiedFluidFSISurface& s1 = (np == 0? m_ss : m_ms);
        FETiedFluidFSISurface& s2 = (np == 0? m_ms : m_ss);
        
        // loop over all primary elements
        for (int i=0; i<s1.Elements(); ++i)
        {
            // get the primary surface element
            FESurfaceElement& se1 = s1.Element(i);
            
            // get the nr of nodes and integration points
            int neln1 = se1.Nodes();
            int nint1 = se1.GaussPoints();
            
            // copy the LM vector; we'll need it later
            s1.UnpackLM(se1, LM1);
            
            // we calculate all the metrics we need before we
            // calculate the nodal forces
            for (int j=0; j<nint1; ++j)
            {
                // get the base vectors
                vec3d g[2];
                s1.CoBaseVectors(se1, j, g);
                
                // jacobians: J = |g0xg1|
                n1[j] = g[0] ^ g[1];
                detJ[j] = n1[j].unit();
                
                // integration weights
                w[j] = se1.GaussWeights()[j];
            }
            
            // loop over all integration points
            // note that we are integrating over the current surface
            for (int j=0; j<nint1; ++j)
            {
                FETiedFluidFSISurface::Data& pt = static_cast<FETiedFluidFSISurface::Data&>(*se1.GetMaterialPoint(j));

                // get the secondary element
                FESurfaceElement* pme = pt.m_pme;
                if (pme)
                {
                    // get the secondary element
                    FESurfaceElement& se2 = *pme;
                    
                    // get the nr of secondary element nodes
                    int neln2 = se2.Nodes();
                    
                    // copy LM vector
                    s2.UnpackLM(se2, LM2);
                    
                    mat3ds sf = GetFluidStress(se2, pt.m_rs);
                    
                    // calculate degrees of freedom
                    int ndof = 3*(neln1 + neln2);
                    
                    // build the LM vector
                    LM.resize(ndof);
                    for (int a=0; a<neln1; ++a)
                    {
                        LM[3*a  ] = LM1[3*a  ];
                        LM[3*a+1] = LM1[3*a+1];
                        LM[3*a+2] = LM1[3*a+2];
                    }
                    
                    for (int b=0; b<neln2; ++b)
                    {
                        LM[3*(b+neln1)  ] = LM2[3*b  ];
                        LM[3*(b+neln1)+1] = LM2[3*b+1];
                        LM[3*(b+neln1)+2] = LM2[3*b+2];
                    }
                    
                    // build the en vector
                    en.resize(neln1+neln2);
                    for (int a=0; a<neln1; ++a) en[a      ] = se1.m_node[a];
                    for (int b=0; b<neln2; ++b) en[b+neln1] = se2.m_node[b];
                    
                    // get primary element shape functions
                    N1 = se1.H(j);
                    
                    // get secondary element shape functions
                    double r = pt.m_rs[0];
                    double s = pt.m_rs[1];
                    se2.shape_fnc(N2, r, s);
                    
                    // gap function
                    vec3d dg = pt.m_dg;
                    
                    // lagrange multiplier
                    vec3d Lm = pt.m_Lmd;
                    
                    // penalty
                    double eps = m_epsn*pt.m_epsn;
                    
                    // contact traction
                    vec3d t1 = Lm + dg*eps - sf*n1[j];
                    pt.m_tr = t1;
                    
                    // calculate the force vector
                    fe.resize(ndof);
                    zero(fe);
                    
                    for (int a=0; a<neln1; ++a) f1[a] = t1*N1[a];
                    for (int b=0; b<neln2; ++b) f2[b] = -t1*N2[b];
                    
                    for (int a = 0; a<neln1; ++a) {
                        fe[3*a  ] += f1[a].x*detJ[j]*w[j];
                        fe[3*a+1] += f1[a].y*detJ[j]*w[j];
                        fe[3*a+2] += f1[a].z*detJ[j]*w[j];
                    }
                    for (int b = 0; b<neln2; ++b) {
                        fe[3*(b+neln1)  ] += f2[b].x*detJ[j]*w[j];
                        fe[3*(b+neln1)+1] += f2[b].y*detJ[j]*w[j];
                        fe[3*(b+neln1)+2] += f2[b].z*detJ[j]*w[j];
                    }

                    for (int a = 0; a < neln1; ++a) s1.m_Fn[se1.m_lnode[a]] += f1[a];
                    for (int b = 0; b < neln2; ++b) s2.m_Fn[se2.m_lnode[b]] += f2[b];

                    // assemble the global residual
                    R.Assemble(en, LM, fe);
                }
            }
        }
    }
}

//-----------------------------------------------------------------------------
void FETiedFluidFSI::StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp)
{
    vector<int> LM1, LM2, LM, en;
    const int MI = FEElement::MAX_INTPOINTS;
    const int MN = FEElement::MAX_NODES;
    double detJ[MI], w[MI], *N1, N2[MI];
    FEElementMatrix ke;
    
    // see how many reformations we've had to do so far
    int nref = GetSolver()->m_nref;
    
    // set higher order stiffness mutliplier
    // NOTE: this algrotihm doesn't really need this
    // but I've added this functionality to compare with the other contact
    // algorithms and to see the effect of the different stiffness contributions
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
    
    double dt = tp.timeIncrement;
    double alpha = tp.alphaf;
    double ksi = tp.gamma / (tp.beta*dt);

    // do single- or two-pass
    int npass = 1;
    for (int np=0; np < npass; ++np)
    {
        // get the primary and secondary surface
        FETiedFluidFSISurface& s1 = (np == 0? m_ss : m_ms);
        FETiedFluidFSISurface& s2 = (np == 0? m_ms : m_ss);
        
        // loop over all primary elements
        for (int i=0; i<s1.Elements(); ++i)
        {
            // get ths primary element
            FESurfaceElement& se1 = s1.Element(i);
            
            // get nr of nodes and integration points
            int neln1 = se1.Nodes();
            int nint1 = se1.GaussPoints();
            
            // copy the LM vector
            s1.UnpackLM(se1, LM1);
            
            // we calculate all the metrics we need before we
            // calculate the nodal forces
            for (int j=0; j<nint1; ++j)
            {
                // get the base vectors
                vec3d g[2];
                s1.CoBaseVectors(se1, j, g);
                
                // jacobians: J = |g0xg1|
                detJ[j] = (g[0] ^ g[1]).norm();
                
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
                    FEElement* e2 = se2.m_elem[0].pe;
                    
                    // get the nr of secondary nodes
                    int neln2 = se2.Nodes();
                    
                    int nint2 = se2.GaussPoints();
                    
                    // copy the LM vector
                    s2.UnpackLM(se2, LM2);
                    
                    int ndpn;    // number of dofs per node
                    int ndof;    // number of dofs in stiffness matrix

                    FEFluidFSI* pfsi = dynamic_cast<FEFluidFSI*>(GetFEModel()->GetMaterial(e2->GetMatID()));
                    assert(pfsi);
                    mat3d Ls2; Ls2.zero();
                    tens4ds Cv2; Cv2.zero();
                    for (int n = 0; n<e2->GaussPoints(); ++n)
                    {
                        FEMaterialPoint& mp = *e2->GetMaterialPoint(n);
                        FEElasticMaterialPoint& ep = *(mp.ExtractData<FEElasticMaterialPoint>());
                        Ls2 += ep.m_L;
                        Cv2 += pfsi->Fluid()->GetViscous()->Tangent_RateOfDeformation(mp);
                    }
                    Ls2 /= e2->GaussPoints();
                    Cv2 /= e2->GaussPoints();

                    // get the fluid stress
                    mat3ds sf = GetFluidStress(se2, pt.m_rs);
                    
                    // calculate degrees of freedom for elastic-on-elastic contact
                    ndpn = 3;
                    ndof = ndpn*(neln1 + neln2);
                    
                    // build the LM vector
                    LM.resize(ndof);
                    
                    for (int a=0; a<neln1; ++a)
                    {
                        LM[3*a  ] = LM1[3*a  ];
                        LM[3*a+1] = LM1[3*a+1];
                        LM[3*a+2] = LM1[3*a+2];
                    }
                    
                    for (int b=0; b<neln2; ++b)
                    {
                        LM[3*(b+neln1)  ] = LM2[3*b  ];
                        LM[3*(b+neln1)+1] = LM2[3*b+1];
                        LM[3*(b+neln1)+2] = LM2[3*b+2];
                    }

                    // build the en vector
                    en.resize(neln1+neln2);
                    for (int a=0; a<neln1; ++a) en[a      ] = se1.m_node[a];
                    for (int b=0; b<neln2; ++b) en[b+neln1] = se2.m_node[b];
                    
                    // primary shape functions
                    N1 = se1.H(j);
                    
                    // secondary shape functions
                    double r = pt.m_rs[0];
                    double s = pt.m_rs[1];
                    se2.shape_fnc(N2, r, s);
                    
                    // get primary normal vector
                    vec3d n1 = pt.m_nu;
                    
                    // gap function
                    vec3d dg = pt.m_dg;
                    
                    // lagrange multiplier
                    vec3d Lm = pt.m_Lmd;
                    
                    // penalty
                    double eps = m_epsn*pt.m_epsn;
                    
                    // contact traction
                    vec3d t1 = Lm + dg*eps - sf*n1;
                    
                    // create the stiffness matrix
                    ke.resize(ndof, ndof); ke.zero();
                    
                    // --- S O L I D - S O L I D   C O N T A C T ---
                    
                    // a. I-term
                    //------------------------------------
                    
                    for (int a=0; a<neln1; ++a) {
                        for (int c=0; c<neln1; ++c)
                        {
                            mat3dd K11(eps*N1[a]*N1[c]*detJ[j]*w[j]*alpha);
                            ke[ndpn*a    ][ndpn*c    ] += K11.xx();
                            ke[ndpn*a + 1][ndpn*c + 1] += K11.yy();
                            ke[ndpn*a + 2][ndpn*c + 2] += K11.zz();
                        }
                        for (int d=0; d<neln2; ++d)
                        {
                            mat3dd K12(-eps*N1[a]*N2[d]*detJ[j]*w[j]*alpha);
                            ke[ndpn*a    ][ndpn*(neln1+d)    ] += K12.xx();
                            ke[ndpn*a + 1][ndpn*(neln1+d) + 1] += K12.yy();
                            ke[ndpn*a + 2][ndpn*(neln1+d) + 2] += K12.zz();
                        }
                    }

                    for (int b=0; b<neln2; ++b) {
                        for (int c=0; c<neln1; ++c)
                        {
                            mat3dd K21(-eps*N2[b]*N1[c]*detJ[j]*w[j]*alpha);
                            ke[ndpn*(neln1+b)    ][ndpn*c    ] += K21.xx();
                            ke[ndpn*(neln1+b) + 1][ndpn*c + 1] += K21.yy();
                            ke[ndpn*(neln1+b) + 2][ndpn*c + 2] += K21.zz();
                        }
                        for (int d=0; d<neln2; ++d)
                        {
                            mat3dd K22(eps*N2[b]*N2[d]*detJ[j]*w[j]*alpha);
                            ke[ndpn*(neln1+b)    ][ndpn*(neln1+d)    ] += K22.xx();
                            ke[ndpn*(neln1+b) + 1][ndpn*(neln1+d) + 1] += K22.yy();
                            ke[ndpn*(neln1+b) + 2][ndpn*(neln1+d) + 2] += K22.zz();
                        }
                    }
 
                    // b. A-term
                    //-------------------------------------
                    
                    double* Gr1 = se1.Gr(j);
                    double* Gs1 = se1.Gs(j);
                    vec3d g1[2];
                    s1.CoBaseVectors(se1, j, g1);
                    
                    vec3d a1[MN];
                    mat3d A1[MN];
                    for (int c=0; c<neln1; ++c) {
                        a1[c] = n1 ^ (g1[1]*Gr1[c] - g1[0]*Gs1[c]);
                        A1[c] = t1 & a1[c];
                    }
                    
                    double Gr2[MN], Gs2[MN];
                    vec3d g2[2];
                    s2.ContraBaseVectors(se2, pt.m_rs.x(), pt.m_rs.y(), g2);
                    se2.shape_deriv(Gr2, Gs2, pt.m_rs.x(), pt.m_rs.y());
                    vec3d gradN2[MN];
                    for (int d=0; d<neln2; ++d) gradN2[d] = g2[0]*Gr2[d] + g2[1]*Gs2[d];
                    
                    if (!m_bsymm)
                    {
                        // non-symmetric
                        for (int a=0; a<neln1; ++a) {
                            for (int c=0; c<neln1; ++c)
                            {
                                mat3d A11 = A1[c]*(N1[a]*w[j]*alpha);
                                ke[ndpn*a    ][ndpn*c    ] += A11(0,0);
                                ke[ndpn*a    ][ndpn*c + 1] += A11(0,1);
                                ke[ndpn*a    ][ndpn*c + 2] += A11(0,2);
                                
                                ke[ndpn*a + 1][ndpn*c    ] += A11(1,0);
                                ke[ndpn*a + 1][ndpn*c + 1] += A11(1,1);
                                ke[ndpn*a + 1][ndpn*c + 2] += A11(1,2);
                                
                                ke[ndpn*a + 2][ndpn*c    ] += A11(2,0);
                                ke[ndpn*a + 2][ndpn*c + 1] += A11(2,1);
                                ke[ndpn*a + 2][ndpn*c + 2] += A11(2,2);
                            }
                            
                            for (int d=0; d<neln2; ++d)
                            {
                                mat3d H12 = -Cv2.dot(mat3dd(alpha*ksi)-Ls2)*(N1[a]*(gradN2[d]*(g1[0]^g1[1])));
                                ke[ndpn*a    ][ndpn*(neln1+d)    ] += H12(0,0);
                                ke[ndpn*a    ][ndpn*(neln1+d) + 1] += H12(0,1);
                                ke[ndpn*a    ][ndpn*(neln1+d) + 2] += H12(0,2);

                                ke[ndpn*a + 1][ndpn*(neln1+d)    ] += H12(1,0);
                                ke[ndpn*a + 1][ndpn*(neln1+d) + 1] += H12(1,1);
                                ke[ndpn*a + 1][ndpn*(neln1+d) + 2] += H12(1,2);

                                ke[ndpn*a + 2][ndpn*(neln1+d)    ] += H12(2,0);
                                ke[ndpn*a + 2][ndpn*(neln1+d) + 1] += H12(2,1);
                                ke[ndpn*a + 2][ndpn*(neln1+d) + 2] += H12(2,2);
                            }
                        }
                        
                        for (int b=0; b<neln2; ++b) {
                            for (int c=0; c<neln1; ++c)
                            {
                                mat3d A21 = A1[c]*(-N2[b]*w[j]*alpha);
                                ke[ndpn*(neln1+b)    ][ndpn*c    ] += A21(0,0);
                                ke[ndpn*(neln1+b)    ][ndpn*c + 1] += A21(0,1);
                                ke[ndpn*(neln1+b)    ][ndpn*c + 2] += A21(0,2);
                                
                                ke[ndpn*(neln1+b) + 1][ndpn*c    ] += A21(1,0);
                                ke[ndpn*(neln1+b) + 1][ndpn*c + 1] += A21(1,1);
                                ke[ndpn*(neln1+b) + 1][ndpn*c + 2] += A21(1,2);
                                
                                ke[ndpn*(neln1+b) + 2][ndpn*c    ] += A21(2,0);
                                ke[ndpn*(neln1+b) + 2][ndpn*c + 1] += A21(2,1);
                                ke[ndpn*(neln1+b) + 2][ndpn*c + 2] += A21(2,2);
                            }
                            
                            for (int d=0; d<neln2; ++d)
                            {
                                mat3d H22 = Cv2.dot(mat3dd(alpha*ksi)-Ls2)*(N2[b]*(gradN2[d]*(g1[0]^g1[1])));
                                ke[ndpn*(neln1+b)     ][ndpn*(neln1+d)    ] += H22(0,0);
                                ke[ndpn*(neln1+b)     ][ndpn*(neln1+d) + 1] += H22(0,1);
                                ke[ndpn*(neln1+b)     ][ndpn*(neln1+d) + 2] += H22(0,2);

                                ke[ndpn*(neln1+b) + 1][ndpn*(neln1+d)    ] += H22(1,0);
                                ke[ndpn*(neln1+b) + 1][ndpn*(neln1+d) + 1] += H22(1,1);
                                ke[ndpn*(neln1+b) + 1][ndpn*(neln1+d) + 2] += H22(1,2);

                                ke[ndpn*(neln1+b) + 2][ndpn*(neln1+d)    ] += H22(2,0);
                                ke[ndpn*(neln1+b) + 2][ndpn*(neln1+d) + 1] += H22(2,1);
                                ke[ndpn*(neln1+b) + 2][ndpn*(neln1+d) + 2] += H22(2,2);
                            }
                        }
                        
                    }
                    else
                    {
                        // symmetric
                        for (int a=0; a<neln1; ++a) {
                            for (int c=0; c<neln1; ++c)
                            {
                                mat3ds A11 = (A1[c]*(N1[a]*w[j]*alpha)).sym();
                                ke[ndpn*a    ][ndpn*c    ] += A11(0,0);
                                ke[ndpn*a    ][ndpn*c + 1] += A11(0,1);
                                ke[ndpn*a    ][ndpn*c + 2] += A11(0,2);
                                
                                ke[ndpn*a + 1][ndpn*c    ] += A11(1,0);
                                ke[ndpn*a + 1][ndpn*c + 1] += A11(1,1);
                                ke[ndpn*a + 1][ndpn*c + 2] += A11(1,2);
                                
                                ke[ndpn*a + 2][ndpn*c    ] += A11(2,0);
                                ke[ndpn*a + 2][ndpn*c + 1] += A11(2,1);
                                ke[ndpn*a + 2][ndpn*c + 2] += A11(2,2);
                            }
                        }
                        
                        for (int b=0; b<neln2; ++b) {
                            for (int c=0; c<neln1; ++c)
                            {
                                mat3ds A21 = (A1[c]*(-N2[b]*w[j]*alpha)).sym();
                                ke[ndpn*(neln1+b)    ][ndpn*c    ] += A21(0,0);
                                ke[ndpn*(neln1+b)    ][ndpn*c + 1] += A21(0,1);
                                ke[ndpn*(neln1+b)    ][ndpn*c + 2] += A21(0,2);
                                
                                ke[ndpn*(neln1+b) + 1][ndpn*c    ] += A21(1,0);
                                ke[ndpn*(neln1+b) + 1][ndpn*c + 1] += A21(1,1);
                                ke[ndpn*(neln1+b) + 1][ndpn*c + 2] += A21(1,2);
                                
                                ke[ndpn*(neln1+b) + 2][ndpn*c    ] += A21(2,0);
                                ke[ndpn*(neln1+b) + 2][ndpn*c + 1] += A21(2,1);
                                ke[ndpn*(neln1+b) + 2][ndpn*c + 2] += A21(2,2);
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
}

//-----------------------------------------------------------------------------
bool FETiedFluidFSI::Augment(int naug, const FETimeInfo& tp)
{
    // make sure we need to augment
    if (m_laugon != FECore::AUGLAG_METHOD) return true;
    
    int i;
    vec3d Ln;
    bool bconv = true;
    
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
            FETiedFluidFSISurface::Data& ds = static_cast<FETiedFluidFSISurface::Data&>(*el.GetMaterialPoint(j));
			normL0 += ds.m_Lmd*ds.m_Lmd;
        }
    }
    for (int i=0; i<NM; ++i)
    {
		FESurfaceElement& el = m_ms.Element(i);
		for (int j=0; j<el.GaussPoints(); ++j)
        {
            FETiedFluidFSISurface::Data& dm = static_cast<FETiedFluidFSISurface::Data&>(*el.GetMaterialPoint(j));
			normL0 += dm.m_Lmd*dm.m_Lmd;
        }
    }
    
    // b. gap component
    // (is calculated during update)
    double maxgap = 0;
    
    // update Lagrange multipliers
    double normL1 = 0, eps;
    for (i=0; i<NS; ++i)
    {
		FESurfaceElement& el = m_ss.Element(i);
		for (int j = 0; j<el.GaussPoints(); ++j)
		{
            FETiedFluidFSISurface::Data& ds = static_cast<FETiedFluidFSISurface::Data&>(*el.GetMaterialPoint(j));

            // update Lagrange multipliers on primary surface
            eps = m_epsn*ds.m_epsn;
            ds.m_Lmd = ds.m_Lmd + ds.m_dg*eps;
            
            normL1 += ds.m_Lmd*ds.m_Lmd;
            
            maxgap = max(maxgap,sqrt(ds.m_dg*ds.m_dg));
        }
    }
    
	for (int i = 0; i<NM; ++i)
	{
		FESurfaceElement& el = m_ms.Element(i);
		for (int j = 0; j<el.GaussPoints(); ++j)
		{
            FETiedFluidFSISurface::Data& dm = static_cast<FETiedFluidFSISurface::Data&>(*el.GetMaterialPoint(j));

            // update Lagrange multipliers on secondary surface
            eps = m_epsn*dm.m_epsn;
            dm.m_Lmd = dm.m_Lmd + dm.m_dg*eps;
            
            normL1 += dm.m_Lmd*dm.m_Lmd;
            
            maxgap = max(maxgap,sqrt(dm.m_dg*dm.m_dg));
        }
    }
    
    // Ideally normP should be evaluated from the fluid pressure at the
    // contact interface (not easily accessible).  The next best thing
    // is to use the contact traction.
    normP = normL1;
    
    // calculate relative norms
    double lnorm = (normL1 != 0 ? fabs((normL1 - normL0) / normL1) : fabs(normL1 - normL0));
    double pnorm = (normP != 0 ? (normDP/normP) : normDP);
    
    // check convergence
    if ((m_gtol > 0) && (maxgap > m_gtol)) bconv = false;
    
    if ((m_atol > 0) && (lnorm > m_atol)) bconv = false;
    if ((m_atol > 0) && (pnorm > m_atol)) bconv = false;
    
    if (naug < m_naugmin ) bconv = false;
    if (naug >= m_naugmax) bconv = true;
    
    feLog(" sliding interface # %d\n", GetID());
    feLog("                        CURRENT        REQUIRED\n");
    feLog("    D multiplier : %15le", lnorm); if (m_atol > 0) feLog("%15le\n", m_atol); else feLog("       ***\n");
    
    feLog("    maximum gap  : %15le", maxgap);
    if (m_gtol > 0) feLog("%15le\n", m_gtol); else feLog("       ***\n");
    
    return bconv;
}

//-----------------------------------------------------------------------------
void FETiedFluidFSI::Serialize(DumpStream &ar)
{
    // store contact data
    FEContactInterface::Serialize(ar);
    
    // store contact surface data
    m_ms.Serialize(ar);
    m_ss.Serialize(ar);
    
    // serialize pointers
    if (ar.IsShallow() == false)
    {
		SerializeElementPointers(m_ss, m_ms, ar);
    }
}
