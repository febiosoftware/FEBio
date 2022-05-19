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
#include "FETiedMultiphasicInterface.h"
#include "FEBiphasic.h"
#include "FEBiphasicSolute.h"
#include "FETriphasic.h"
#include "FEMultiphasic.h"
#include "FECore/FEModel.h"
#include "FECore/FEAnalysis.h"
#include "FECore/FENormalProjection.h"
#include <FECore/FELinearSystem.h>
#include "FECore/log.h"

//-----------------------------------------------------------------------------
// Define sliding interface parameters
BEGIN_FECORE_CLASS(FETiedMultiphasicInterface, FEContactInterface)
	ADD_PARAMETER(m_atol     , "tolerance"          );
	ADD_PARAMETER(m_gtol     , "gaptol"             );
	ADD_PARAMETER(m_ptol     , "ptol"               );
    ADD_PARAMETER(m_ctol     , "ctol"               );
	ADD_PARAMETER(m_epsn     , "penalty"            );
	ADD_PARAMETER(m_bautopen , "auto_penalty"       );
    ADD_PARAMETER(m_bupdtpen , "update_penalty"     );
	ADD_PARAMETER(m_btwo_pass, "two_pass"           );
	ADD_PARAMETER(m_knmult   , "knmult"             );
	ADD_PARAMETER(m_stol     , "search_tol"         );
	ADD_PARAMETER(m_epsp     , "pressure_penalty"   );
	ADD_PARAMETER(m_epsc     , "concentration_penalty");
	ADD_PARAMETER(m_bsymm    , "symmetric_stiffness");
	ADD_PARAMETER(m_srad     , "search_radius"      );
	ADD_PARAMETER(m_naugmin  , "minaug"             );
	ADD_PARAMETER(m_naugmax  , "maxaug"             );
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FETiedMultiphasicSurface::Data::Data()
{
    m_Gap = vec3d(0,0,0);
    m_dg = vec3d(0,0,0);
    m_nu = vec3d(0,0,0);
    m_rs = vec2d(0,0);
    m_Lmd = vec3d(0,0,0);
    m_Lmp  = 0.0;
    m_epsn = 1.0;
    m_epsp = 1.0;
    m_pg   = 0.0;
}

void FETiedMultiphasicSurface::Data::Serialize(DumpStream& ar)
{
	FEBiphasicContactPoint::Serialize(ar);
	ar & m_Gap;
	ar & m_dg;
	ar & m_nu;
	ar & m_rs;
	ar & m_Lmd;
	ar & m_epsn;
	ar & m_epsp;
    ar & m_Lmc;
	ar & m_epsc;
    ar & m_cg;
}

//-----------------------------------------------------------------------------
// FETiedMultiphasicSurface
//-----------------------------------------------------------------------------

FETiedMultiphasicSurface::FETiedMultiphasicSurface(FEModel* pfem) : FEBiphasicContactSurface(pfem)
{
    m_bporo = m_bsolu = false;
    m_dofC = -1;
}

//-----------------------------------------------------------------------------
//! create material point data
FEMaterialPoint* FETiedMultiphasicSurface::CreateMaterialPoint()
{
	return new FETiedMultiphasicSurface::Data;
}

//-----------------------------------------------------------------------------
void FETiedMultiphasicSurface::UnpackLM(FEElement& el, vector<int>& lm)
{
    // get nodal DOFS
    DOFS& dofs = GetFEModel()->GetDOFS();
    int MAX_CDOFS = dofs.GetVariableSize("concentration");
    
    int N = el.Nodes();
    lm.resize(N*(4+MAX_CDOFS));
    
    // pack the equation numbers
    for (int i=0; i<N; ++i)
    {
        int n = el.m_node[i];
        
        FENode& node = m_pMesh->Node(n);
        vector<int>& id = node.m_ID;
        
        // first the displacement dofs
        lm[3*i  ] = id[m_dofX];
        lm[3*i+1] = id[m_dofY];
        lm[3*i+2] = id[m_dofZ];
        
        // now the pressure dofs
        lm[3*N+i] = id[m_dofP];
        
        // concentration dofs
        for (int k=0; k<MAX_CDOFS; ++k)
            lm[(4 + k)*N + i] = id[m_dofC+k];
    }
}

//-----------------------------------------------------------------------------
bool FETiedMultiphasicSurface::Init()
{
    // initialize surface data first
    if (FEBiphasicContactSurface::Init() == false) return false;
    
    // store concentration index
    DOFS& dofs = GetFEModel()->GetDOFS();
    m_dofC = dofs.GetDOF("concentration", 0);
    
    // allocate node normals
    m_nn.assign(Nodes(), vec3d(0,0,0));
    
    // determine solutes for this surface using the first surface element
    // TODO: Check that all elements use the same set of solutes as the first element
    int nsol = 0;
    if (Elements()) {
        FESurfaceElement& se = Element(0);
        // get the element this surface element belongs to
        FEElement* pe = se.m_elem[0];
        if (pe)
        {
            // get the material
            FEMaterial* pm = m_pfem->GetMaterial(pe->GetMatID());
            
            // check type of element
            FEBiphasic* pb = dynamic_cast<FEBiphasic*> (pm);
            FEBiphasicSolute* pbs = dynamic_cast<FEBiphasicSolute*> (pm);
            FETriphasic* pt = dynamic_cast<FETriphasic*>(pm);
            FEMultiphasic* pmp = dynamic_cast<FEMultiphasic*> (pm);
            if (pb)
            {
                m_bporo = true;
                nsol = 0;
            }
            else if (pbs)
            {
                m_bporo = m_bsolu = true;
                nsol = 1;
                m_sid.assign(nsol, pbs->GetSolute()->GetSoluteID() - 1);
            }
            else if (pt)
            {
                m_bporo = m_bsolu = true;
                nsol = 2;
                m_sid.resize(nsol);
                m_sid[0] = pt->m_pSolute[0]->GetSoluteID() - 1;
                m_sid[1] = pt->m_pSolute[1]->GetSoluteID() - 1;
            }
            else if (pmp)
            {
                m_bporo = m_bsolu = true;
                nsol = pmp->Solutes();
                m_sid.resize(nsol);
                for (int isol=0; isol<nsol; ++isol) {
                    m_sid[isol] = pmp->GetSolute(isol)->GetSoluteID() - 1;
                }
            }
        }
    }
    
    // allocate data structures
    int NE = Elements();
    m_poro.resize(NE,false);
    for (int i=0; i<NE; ++i)
    {
        FESurfaceElement& el = Element(i);
        // get the element this surface element belongs to
        FEElement* pe = el.m_elem[0];
        if (pe)
        {
            // get the material
            FEMaterial* pm = m_pfem->GetMaterial(pe->GetMatID());
            
            // see if this is a poro-elastic element
            FEBiphasic* bp = dynamic_cast<FEBiphasic*> (pm);
            FEBiphasicSolute* bs = dynamic_cast<FEBiphasicSolute*> (pm);
            FETriphasic* tp = dynamic_cast<FETriphasic*> (pm);
            FEMultiphasic* mp = dynamic_cast<FEMultiphasic*> (pm);
            if (bp || bs || tp || mp) {
                m_poro[i] = true;
                m_bporo = true;
            }
            if (bs || tp || mp) {
                m_bsolu = true;
            }
        }
        int nint = el.GaussPoints();
        if (nsol) {
            for (int j=0; j<nint; ++j) {
				Data& data = static_cast<Data&>(*el.GetMaterialPoint(j));
                data.m_Lmc.resize(nsol);
                data.m_epsc.resize(nsol);
                data.m_cg.resize(nsol);
            }
        }
    }
    
    return true;
}

//-----------------------------------------------------------------------------
//! This function calculates the node normal. Due to the piecewise continuity
//! of the surface elements this normal is not uniquely defined so in order to
//! obtain a unique normal the normal is averaged for each node over all the
//! element normals at the node

void FETiedMultiphasicSurface::UpdateNodeNormals()
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
void FETiedMultiphasicSurface::Serialize(DumpStream& ar)
{
	FEBiphasicContactSurface::Serialize(ar);
	ar & m_bporo;
	ar & m_bsolu;
	ar & m_poro;
	ar & m_nn;
	ar & m_sid;
}

//-----------------------------------------------------------------------------
// FETiedMultiphasicInterface
//-----------------------------------------------------------------------------

FETiedMultiphasicInterface::FETiedMultiphasicInterface(FEModel* pfem) : FEContactInterface(pfem), m_ss(pfem), m_ms(pfem)
{
    static int count = 1;
    SetID(count++);
    
    // initial values
    m_knmult = 1;
    m_atol = 0.1;
    m_epsn = 1;
    m_epsp = 1;
    m_epsc = 1;
    m_btwo_pass = false;
    m_stol = 0.01;
    m_bsymm = true;
    m_srad = 1.0;
    m_gtol = 0;
    m_ptol = 0;
    m_ctol = 0;
    m_bautopen = false;
    m_bupdtpen = false;
    
    m_naugmin = 0;
    m_naugmax = 10;
    
    m_dofP = pfem->GetDOFIndex("p");
    m_dofC = pfem->GetDOFIndex("concentration", 0);
    
    // set parents
    m_ss.SetContactInterface(this);
    m_ms.SetContactInterface(this);

    m_ss.SetSibling(&m_ms);
    m_ms.SetSibling(&m_ss);
}

//-----------------------------------------------------------------------------

FETiedMultiphasicInterface::~FETiedMultiphasicInterface()
{
}

//-----------------------------------------------------------------------------
bool FETiedMultiphasicInterface::Init()
{
    m_Rgas = GetFEModel()->GetGlobalConstant("R");
    m_Tabs = GetFEModel()->GetGlobalConstant("T");
    
    // initialize surface data
    if (m_ss.Init() == false) return false;
    if (m_ms.Init() == false) return false;
    
    // determine which solutes are common to both contact surfaces
    m_sid.clear(); m_ssl.clear(); m_msl.clear(); m_sz.clear();
    for (int is=0; is<m_ss.m_sid.size(); ++is) {
        for (int im=0; im<m_ms.m_sid.size(); ++im) {
            if (m_ms.m_sid[im] == m_ss.m_sid[is]) {
                m_sid.push_back(m_ss.m_sid[is]);
                m_ssl.push_back(is);
                m_msl.push_back(im);
                FESoluteData* sd = FindSoluteData(m_ss.m_sid[is]);
                m_sz.push_back(sd->m_z);
            }
        }
    }
    
    return true;
}

//-----------------------------------------------------------------------------
//! build the matrix profile for use in the stiffness matrix
void FETiedMultiphasicInterface::BuildMatrixProfile(FEGlobalMatrix& K)
{
    FEModel& fem = *GetFEModel();
    FEMesh& mesh = fem.GetMesh();
    
    // get the DOFS
    const int dof_X = fem.GetDOFIndex("x");
    const int dof_Y = fem.GetDOFIndex("y");
    const int dof_Z = fem.GetDOFIndex("z");
    const int dof_P = fem.GetDOFIndex("p");
    const int dof_C = fem.GetDOFIndex("concentration", 0);
    
    int nsol = (int)m_sid.size();
    int ndpn = 7 + nsol;
    
    vector<int> lm(ndpn*FEElement::MAX_NODES*2);
    
    int npass = (m_btwo_pass?2:1);
    for (int np=0; np<npass; ++np)
    {
        FETiedMultiphasicSurface& ss = (np == 0? m_ss : m_ms);
        
        int k, l;
        for (int j=0; j<ss.Elements(); ++j)
        {
            FESurfaceElement& se = ss.Element(j);
            int nint = se.GaussPoints();
            int* sn = &se.m_node[0];
            for (k=0; k<nint; ++k)
            {
				FETiedMultiphasicSurface::Data& data = static_cast<FETiedMultiphasicSurface::Data&>(*se.GetMaterialPoint(k));
                FESurfaceElement* pe = data.m_pme;
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
                        lm[ndpn*l  ] = id[dof_X];
                        lm[ndpn*l+1] = id[dof_Y];
                        lm[ndpn*l+2] = id[dof_Z];
                        lm[ndpn*l+3] = id[dof_P];
                        for (int m=0; m<nsol; ++m) {
                            lm[ndpn*l+4+m] = id[dof_C + m_sid[m]];
                        }
                    }
                    
                    for (l=0; l<nmeln; ++l)
                    {
                        vector<int>& id = mesh.Node(mn[l]).m_ID;
                        lm[ndpn*(l+nseln)  ] = id[dof_X];
                        lm[ndpn*(l+nseln)+1] = id[dof_Y];
                        lm[ndpn*(l+nseln)+2] = id[dof_Z];
                        lm[ndpn*(l+nseln)+3] = id[dof_P];
                        for (int m=0; m<nsol; ++m) {
                            lm[ndpn*(l+nseln)+4+m] = id[dof_C + m_sid[m]];
                        }
                    }
                    
                    K.build_add(lm);
                }
            }
        }
    }
    
}

//-----------------------------------------------------------------------------
void FETiedMultiphasicInterface::UpdateAutoPenalty()
{
    // calculate the penalty
    if (m_bautopen)
    {
        CalcAutoPenalty(m_ss);
        CalcAutoPenalty(m_ms);
        if (m_ss.m_bporo) CalcAutoPressurePenalty(m_ss);
        for (int is=0; is<m_ssl.size(); ++is)
            CalcAutoConcentrationPenalty(m_ss, m_ssl[is]);
        if (m_ms.m_bporo) CalcAutoPressurePenalty(m_ms);
        for (int im=0; im<m_msl.size(); ++im)
            CalcAutoConcentrationPenalty(m_ms, m_msl[im]);
    }
}

//-----------------------------------------------------------------------------
void FETiedMultiphasicInterface::Activate()
{
    // don't forget to call the base class
    FEContactInterface::Activate();
    
    UpdateAutoPenalty();
    
    // project the surfaces onto each other
    // this will evaluate the gap functions in the reference configuration
    InitialProjection(m_ss, m_ms);
    if (m_btwo_pass) InitialProjection(m_ms, m_ss);
}

//-----------------------------------------------------------------------------
void FETiedMultiphasicInterface::CalcAutoPenalty(FETiedMultiphasicSurface& s)
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
			FETiedMultiphasicSurface::Data& pt = static_cast<FETiedMultiphasicSurface::Data&>(*el.GetMaterialPoint(j));
			pt.m_epsn = eps;
        }
    }
}

//-----------------------------------------------------------------------------
//! This function calculates a contact penalty parameter based on the
//! material and geometrical properties of the primary and secondary surfaces
//!
double FETiedMultiphasicInterface::AutoPenalty(FESurfaceElement& el, FESurface &s)
{
    // get the mesh
    FEMesh& m = GetFEModel()->GetMesh();
    
    // get the element this surface element belongs to
    FEElement* pe = el.m_elem[0];
    if (pe == 0) return 0.0;

    tens4ds S;
    // get a material point
    FEMaterialPoint& mp = *pe->GetMaterialPoint(0);
        
    // extract the material
    FEMaterial* pme = GetFEModel()->GetMaterial(pe->GetMatID());
	if (pme == 0) return 0.0;
        
	// get the tangent (stiffness)
	if (dynamic_cast<FEMultiphasic*>(pme)) {
		FEMultiphasic* pmm = dynamic_cast<FEMultiphasic*>(pme);
		S = pmm->Tangent(mp);
	}
	else if (dynamic_cast<FETriphasic*>(pme)) {
		FETriphasic* pmt = dynamic_cast<FETriphasic*>(pme);
		S = pmt->Tangent(mp);
	}
	else if (dynamic_cast<FEBiphasicSolute*>(pme)) {
		FEBiphasicSolute* pms = dynamic_cast<FEBiphasicSolute*>(pme);
		S = pms->Tangent(mp);
	}
	else if (dynamic_cast<FEBiphasic*>(pme)) {
		FEBiphasic* pmb = dynamic_cast<FEBiphasic*>(pme);
		S = pmb->Tangent(mp);
	}
	else if (dynamic_cast<FEElasticMaterial*>(pme)) {
		FEElasticMaterial* pm = dynamic_cast<FEElasticMaterial*>(pme);
		S = pm->Tangent(mp);
	}
	// get the inverse (compliance) at this point
	tens4ds C = S.inverse();
            
	// evaluate element surface normal at parametric center
	vec3d t[2];
	s.CoBaseVectors0(el, 0, 0, t);
	vec3d n = t[0] ^ t[1];
	n.unit();
            
	// evaluate normal component of the compliance matrix
	// (equivalent to inverse of Young's modulus along n)
	double eps = 1./(n*(vdotTdotv(n, C, n)*n));
    
	// get the area of the surface element
	double A = s.FaceArea(el);

	// get the volume of the volume element
	double V = m.ElementVolume(*pe);

    return eps*A/V;
}

//-----------------------------------------------------------------------------
void FETiedMultiphasicInterface::CalcAutoPressurePenalty(FETiedMultiphasicSurface& s)
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
			FETiedMultiphasicSurface::Data& pt = static_cast<FETiedMultiphasicSurface::Data&>(*el.GetMaterialPoint(j));
			pt.m_epsp = eps;
        }
    }
}

//-----------------------------------------------------------------------------

double FETiedMultiphasicInterface::AutoPressurePenalty(FESurfaceElement& el, FETiedMultiphasicSurface& s)
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
        
    // get a material point
    FEMaterialPoint& mp = *pe->GetMaterialPoint(0);
        
    mat3ds K;
        
    // check type of element
    FEBiphasic* pb = dynamic_cast<FEBiphasic*> (pm);
    FEBiphasicSolute* pbs = dynamic_cast<FEBiphasicSolute*> (pm);
    FETriphasic* ptp = dynamic_cast<FETriphasic*> (pm);
    FEMultiphasic* pmp = dynamic_cast<FEMultiphasic*> (pm);
    if (pb) {
        double k[3][3];
        pb->Permeability(k, mp);
        K = mat3ds(k[0][0], k[1][1], k[2][2], k[0][1], k[1][2], k[0][2]);
    }
    else if (ptp)
        K = ptp->GetPermeability()->Permeability(mp);
    else if (pbs)
        K = pbs->GetPermeability()->Permeability(mp);
    else if (pmp)
        K = pmp->GetPermeability()->Permeability(mp);
        
    double eps = n*(K*n);
    
	// get the area of the surface element
	double A = s.FaceArea(el);

	// get the volume of the volume element
	double V = m.ElementVolume(*pe);

    return eps*A/V;
}

//-----------------------------------------------------------------------------
void FETiedMultiphasicInterface::CalcAutoConcentrationPenalty(FETiedMultiphasicSurface& s, const int isol)
{
    // loop over all surface elements
    int ni = 0;
    for (int i=0; i<s.Elements(); ++i)
    {
        // get the surface element
        FESurfaceElement& el = s.Element(i);
        
        // calculate a penalty
        double eps = AutoConcentrationPenalty(el, s, isol);
        
        // assign to integation points of surface element
        int nint = el.GaussPoints();
        for (int j=0; j<nint; ++j, ++ni)
        {
			FETiedMultiphasicSurface::Data& pt = static_cast<FETiedMultiphasicSurface::Data&>(*el.GetMaterialPoint(j));
			pt.m_epsc[isol] = eps;
        }
    }
}

//-----------------------------------------------------------------------------
double FETiedMultiphasicInterface::AutoConcentrationPenalty(FESurfaceElement& el, FETiedMultiphasicSurface& s, const int isol)
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
        
    // get a material point
    FEMaterialPoint& mp = *pe->GetMaterialPoint(0);
        
    mat3ds D;
        
    // see if this is a biphasic-solute or multiphasic element
    FEBiphasicSolute* pbs = dynamic_cast<FEBiphasicSolute*> (pm);
    FETriphasic* ptp = dynamic_cast<FETriphasic*> (pm);
    FEMultiphasic* pmp = dynamic_cast<FEMultiphasic*> (pm);
    if (pbs)
    {
        D = pbs->GetSolute()->m_pDiff->Diffusivity(mp)
        *(pbs->Porosity(mp)*pbs->GetSolute()->m_pSolub->Solubility(mp));
    }
    else if (ptp)
    {
        D = ptp->GetSolute(isol)->m_pDiff->Diffusivity(mp)
        *(ptp->Porosity(mp)*ptp->GetSolute(isol)->m_pSolub->Solubility(mp));
    }
    else if (pmp)
    {
        D = pmp->GetSolute(isol)->m_pDiff->Diffusivity(mp)
        *(pmp->Porosity(mp)*pmp->GetSolute(isol)->m_pSolub->Solubility(mp));
    }
        
    // evaluate normal component of diffusivity
	double eps = n*(D*n);
    
	// get the area of the surface element
	double A = s.FaceArea(el);

	// get the volume of the volume element
	double V = m.ElementVolume(*pe);

    return eps*A/V;
}

//-----------------------------------------------------------------------------
// Perform initial projection between tied surfaces in reference configuration
void FETiedMultiphasicInterface::InitialProjection(FETiedMultiphasicSurface& ss, FETiedMultiphasicSurface& ms)
{
    FESurfaceElement* pme;
    vec3d r, nu;
    double rs[2];
    
    // initialize projection data
    FENormalProjection np(ms);
    np.SetTolerance(m_stol);
    np.SetSearchRadius(m_srad);
    np.Init();
    
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
            
			FETiedMultiphasicSurface::Data& pt = static_cast<FETiedMultiphasicSurface::Data&>(*el.GetMaterialPoint(j));
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
            }
            else
            {
                // the node is not in contact
                pt.m_Gap = vec3d(0,0,0);
            }
        }
    }
}

//-----------------------------------------------------------------------------
// Evaluate gap functions for position and fluid pressure
void FETiedMultiphasicInterface::ProjectSurface(FETiedMultiphasicSurface& ss, FETiedMultiphasicSurface& ms)
{
    FEMesh& mesh = GetFEModel()->GetMesh();
    FESurfaceElement* pme;
    vec3d r;
    double rs[2];
    
    const int MN = FEElement::MAX_NODES;
    double ps[MN], p1;
    int nsol = (int)m_sid.size();
    vector< vector<double> > cs(nsol, vector<double>(MN));
    vector<double> c1(nsol);
    
    // loop over all integration points
    for (int i=0; i<ss.Elements(); ++i)
    {
        FESurfaceElement& el = ss.Element(i);
        bool sporo = ss.m_poro[i];
        
        int ne = el.Nodes();
        int nint = el.GaussPoints();
        
        // get the nodal pressures
        p1 = 0;
        if (sporo)
        {
            for (int j=0; j<ne; ++j) ps[j] = mesh.Node(el.m_node[j]).get(m_dofP);
        }
        
        // get the nodal concentrations
        for (int isol=0; isol<nsol; ++isol) {
            for (int j=0; j<ne; ++j) cs[isol][j] = mesh.Node(el.m_node[j]).get(m_dofC + m_sid[isol]);
        }
        
        for (int j=0; j<nint; ++j)
        {
			FETiedMultiphasicSurface::Data& pt = static_cast<FETiedMultiphasicSurface::Data&>(*el.GetMaterialPoint(j));

            // calculate the global position of the integration point
            r = ss.Local2Global(el, j);
            
            // get the pressure at the integration point
            if (sporo) p1 = el.eval(ps, j);
            
            // get the concentration at the integration point
            for (int isol=0; isol<nsol; ++isol) c1[isol] = el.eval(&cs[isol][0], j);
            
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
                
                // calculate the pressure gap function
                bool mporo = ms.m_poro[pme->m_lid];
                if (sporo && mporo) {
                    double pm[FEElement::MAX_NODES];
                    for (int k=0; k<pme->Nodes(); ++k) pm[k] = mesh.Node(pme->m_node[k]).get(m_dofP);
                    double p2 = pme->eval(pm, pt.m_rs[0], pt.m_rs[1]);
                    pt.m_pg = p1 - p2;
                }
                // calculate the concentration gap functions
                for (int isol=0; isol<nsol; ++isol) {
                    int sid = m_sid[isol];
                    double cm[MN];
                    for (int k=0; k<pme->Nodes(); ++k) cm[k] = mesh.Node(pme->m_node[k]).get(m_dofC + sid);
                    double c2 = pme->eval(cm, rs[0], rs[1]);
                    pt.m_cg[m_ssl[isol]] = c1[isol] - c2;
                }
            }
            else
            {
                // the node is not tied
                pt.m_dg = vec3d(0,0,0);
                if (sporo) {
                    pt.m_pg = 0;
                    pt.m_Lmp = 0;
                }
                for (int isol=0; isol<nsol; ++isol) {
                    pt.m_Lmc[m_ssl[isol]] = 0;
                    pt.m_cg[m_ssl[isol]] = 0;
                }
            }
        }
    }
}

//-----------------------------------------------------------------------------

void FETiedMultiphasicInterface::Update()
{
    // project the surfaces onto each other
    // this will update the gap functions as well
    ProjectSurface(m_ss, m_ms);
    if (m_btwo_pass) ProjectSurface(m_ms, m_ss);
    
}

//-----------------------------------------------------------------------------
void FETiedMultiphasicInterface::LoadVector(FEGlobalVector& R, const FETimeInfo& tp)
{
    int i, j, k;
    vector<int> sLM, mLM, LM, en;
    vector<double> fe;
    const int MN = FEElement::MAX_NODES;
    double detJ[MN], w[MN], *Hs, Hm[MN];
    double N[MN*10];
    int nsol = (int)m_sid.size();
    vec3d tn[MN];
    double wn[MN];
    vector< vector<double> > jn(nsol,vector<double>(MN));
    
    FEModel& fem = *GetFEModel();
    FEAnalysis* pstep = fem.GetCurrentStep();
    FESolver* psolver = pstep->GetFESolver();

    double dt = fem.GetTime().timeIncrement;
    
    // Update auto-penalty if requested
    if (m_bupdtpen && psolver->m_niter == 0) UpdateAutoPenalty();
    
    // loop over the nr of passes
    int npass = (m_btwo_pass?2:1);
    for (int np=0; np<npass; ++np)
    {
        // get primary and seconary surface
        FETiedMultiphasicSurface& ss = (np == 0? m_ss : m_ms);
        FETiedMultiphasicSurface& ms = (np == 0? m_ms : m_ss);
        vector<int>& sl = (np == 0? m_ssl : m_msl);
        
        // loop over all primary surface elements
        for (i=0; i<ss.Elements(); ++i)
        {
            // get the surface element
            FESurfaceElement& se = ss.Element(i);
            
            bool sporo = ss.m_bporo;
            
            // get the nr of nodes and integration points
            int nseln = se.Nodes();
            int nint = se.GaussPoints();
            
            // copy the LM vector; we'll need it later
            ss.UnpackLM(se, sLM);
            
            // we calculate all the metrics we need before we
            // calculate the nodal forces
            for (j=0; j<nint; ++j)
            {
                // get the base vectors
                vec3d g[2];
                ss.CoBaseVectors(se, j, g);
                
                // jacobians: J = |g0xg1|
                detJ[j] = (g[0] ^ g[1]).norm();
                
                // integration weights
                w[j] = se.GaussWeights()[j];
                
				FETiedMultiphasicSurface::Data& pt = static_cast<FETiedMultiphasicSurface::Data&>(*se.GetMaterialPoint(j));

                // contact traction
                double eps = m_epsn*pt.m_epsn;      // penalty
                tn[j] = pt.m_Lmd + pt.m_dg*eps;    // contact traction
                
                // normal fluid flux
                double epsp = m_epsp*pt.m_epsp;
                wn[j] = pt.m_Lmp + epsp*pt.m_pg;
                
                // normal solute flux
                for (int isol=0; isol<nsol; ++isol)
                {
                    int l = sl[isol];
                    double epsc = m_epsc*pt.m_epsc[l];
                    jn[isol][j] = pt.m_Lmc[l] + epsc*pt.m_cg[l];
                }
            }
            
            // loop over all integration points
            // note that we are integrating over the current surface
            for (j=0; j<nint; ++j)
            {
				FETiedMultiphasicSurface::Data& pt = static_cast<FETiedMultiphasicSurface::Data&>(*se.GetMaterialPoint(j));
				// get the secondary surface element
                FESurfaceElement* pme = pt.m_pme;
                if (pme)
                {
                    // get the secondary surface element
                    FESurfaceElement& me = *pme;
                    
                    bool mporo = ms.m_bporo;
                    
                    // get the nr of secondary surface element nodes
                    int nmeln = me.Nodes();
                    
                    // copy LM vector
                    ms.UnpackLM(me, mLM);
                    
                    // calculate degrees of freedom
                    int ndof = 3*(nseln + nmeln);
                    
                    // build the LM vector
                    LM.resize(ndof);
                    for (k=0; k<nseln; ++k)
                    {
                        LM[3*k  ] = sLM[3*k  ];
                        LM[3*k+1] = sLM[3*k+1];
                        LM[3*k+2] = sLM[3*k+2];
                    }
                    
                    for (k=0; k<nmeln; ++k)
                    {
                        LM[3*(k+nseln)  ] = mLM[3*k  ];
                        LM[3*(k+nseln)+1] = mLM[3*k+1];
                        LM[3*(k+nseln)+2] = mLM[3*k+2];
                    }
                    
                    // build the en vector
                    en.resize(nseln+nmeln);
                    for (k=0; k<nseln; ++k) en[k      ] = se.m_node[k];
                    for (k=0; k<nmeln; ++k) en[k+nseln] = me.m_node[k];
                    
                    // get element shape functions
                    Hs = se.H(j);
                    
                    // get secondary surface element shape functions
                    double r = pt.m_rs[0];
                    double s = pt.m_rs[1];
                    me.shape_fnc(Hm, r, s);
                    
                    // calculate the force vector
                    fe.resize(ndof);
                    zero(fe);
                    
                    for (k=0; k<nseln; ++k)
                    {
                        N[3*k  ] = Hs[k]*tn[j].x;
                        N[3*k+1] = Hs[k]*tn[j].y;
                        N[3*k+2] = Hs[k]*tn[j].z;
                    }
                    
                    for (k=0; k<nmeln; ++k)
                    {
                        N[3*(k+nseln)  ] = -Hm[k]*tn[j].x;
                        N[3*(k+nseln)+1] = -Hm[k]*tn[j].y;
                        N[3*(k+nseln)+2] = -Hm[k]*tn[j].z;
                    }
                    
                    for (k=0; k<ndof; ++k) fe[k] += N[k]*detJ[j]*w[j];
                    
                    // assemble the global residual
                    R.Assemble(en, LM, fe);
                    
                    // do the biphasic stuff
                    if (sporo && mporo)
                    {
                        // calculate nr of pressure dofs
                        int ndof = nseln + nmeln;
                        
                        // fill the LM
                        LM.resize(ndof);
                        for (k=0; k<nseln; ++k) LM[k        ] = sLM[3*nseln+k];
                        for (k=0; k<nmeln; ++k) LM[k + nseln] = mLM[3*nmeln+k];
                        
                        // fill the force array
                        fe.resize(ndof);
                        zero(fe);
                        for (k=0; k<nseln; ++k) N[k      ] =  Hs[k];
                        for (k=0; k<nmeln; ++k) N[k+nseln] = -Hm[k];
                        
                        for (k=0; k<ndof; ++k) fe[k] += dt*N[k]*wn[j]*detJ[j]*w[j];
                        
                        // assemble residual
                        R.Assemble(en, LM, fe);
                    }
                    for (int isol=0; isol<nsol; ++isol)
                    {
                        int sid = m_sid[isol];
                        
                        // calculate nr of concentration dofs
                        int ndof = nseln + nmeln;
                        
                        // fill the LM
                        LM.resize(ndof);
                        for (k=0; k<nseln; ++k) LM[k        ] = sLM[(4+sid)*nseln+k];
                        for (k=0; k<nmeln; ++k) LM[k + nseln] = mLM[(4+sid)*nmeln+k];
                        
                        // fill the force array
                        fe.resize(ndof);
                        zero(fe);
                        for (k=0; k<nseln; ++k) N[k      ] =  Hs[k];
                        for (k=0; k<nmeln; ++k) N[k+nseln] = -Hm[k];
                        
                        for (k=0; k<ndof; ++k) fe[k] += dt*N[k]*jn[isol][j]*detJ[j]*w[j];
                        
                        // assemble residual
                        R.Assemble(en, LM, fe);
                    }
                }
            }
        }
    }
}

//-----------------------------------------------------------------------------
void FETiedMultiphasicInterface::StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp)
{
    int i, j, k, l;
    vector<int> sLM, mLM, LM, en;
    const int MN = FEElement::MAX_NODES;
    double detJ[MN], w[MN], *Hs, Hm[MN];
    FEElementMatrix ke;
    int nsol = (int)m_sid.size();
    vec3d tn[MN];
    double wn[MN];
    vector< vector<double> > jn(nsol,vector<double>(MN));
    
    FEModel& fem = *GetFEModel();
    
    // see how many reformations we've had to do so far
    int nref = LS.GetSolver()->m_nref;
    
    // set higher order stiffness mutliplier
    // NOTE: this algorithm doesn't really need this
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
    
    // do single- or two-pass
    int npass = (m_btwo_pass?2:1);
    for (int np=0; np < npass; ++np)
    {
        // get the primary and secondary surface
        FETiedMultiphasicSurface& ss = (np == 0? m_ss : m_ms);
        FETiedMultiphasicSurface& ms = (np == 0? m_ms : m_ss);
        vector<int>& sl = (np == 0? m_ssl : m_msl);
        
        // loop over all primary surface elements
        for (i=0; i<ss.Elements(); ++i)
        {
            // get the next element
            FESurfaceElement& se = ss.Element(i);
            
            bool sporo = ss.m_bporo;
            
            // get nr of nodes and integration points
            int nseln = se.Nodes();
            int nint = se.GaussPoints();
            
            double pn[MN];
            vector< vector<double> >cn(nsol,vector<double>(MN));
            for (j=0; j<nseln; ++j)
            {
                pn[j] = ss.GetMesh()->Node(se.m_node[j]).get(m_dofP);
                for (int isol=0; isol<nsol; ++isol) {
                    cn[isol][j] = ss.GetMesh()->Node(se.m_node[j]).get(m_dofC + m_sid[isol]);
                }
            }
            
            // copy the LM vector
            ss.UnpackLM(se, sLM);
            
            // we calculate all the metrics we need before we
            // calculate the nodal forces
            for (j=0; j<nint; ++j)
            {
                // get the base vectors
                vec3d g[2];
                ss.CoBaseVectors(se, j, g);
                
                // jacobians: J = |g0xg1|
                detJ[j] = (g[0] ^ g[1]).norm();
                
                // integration weights
                w[j] = se.GaussWeights()[j];
                
				FETiedMultiphasicSurface::Data& pd = static_cast<FETiedMultiphasicSurface::Data&>(*se.GetMaterialPoint(j));

                // contact traction
                double eps = m_epsn*pd.m_epsn;      // penalty
                tn[j] = pd.m_Lmd + pd.m_dg*eps;     // contact traction
                
                // normal fluid flux
                double epsp = m_epsp*pd.m_epsp;
                wn[j] = pd.m_Lmp + epsp*pd.m_pg;
                
                // normal solute flux
                for (int isol=0; isol<nsol; ++isol)
                {
                    int l = sl[isol];
                    double epsc = m_epsc*pd.m_epsc[l];
                    jn[isol][j] = pd.m_Lmc[l] + epsc*pd.m_cg[l];
                }
                
                // contravariant basis vectors of primary surface
                vec3d Gs[2];
                ss.ContraBaseVectors(se, j, Gs);
            }
            
            // loop over all integration points
            for (j=0; j<nint; ++j)
            {
				FETiedMultiphasicSurface::Data& pt = static_cast<FETiedMultiphasicSurface::Data&>(*se.GetMaterialPoint(j));

				// get the secondary surface element
                FESurfaceElement* pme = pt.m_pme;
                if (pme)
                {
                    FESurfaceElement& me = *pme;
                    
                    bool mporo = ms.m_bporo;
                    
                    // get the nr of secondary surface nodes
                    int nmeln = me.Nodes();
                    
                    // nodal data
                    double pm[MN];
                    vector< vector<double> > cm(nsol,vector<double>(MN));
                    for (k=0; k<nmeln; ++k)
                    {
                        pm[k] = ms.GetMesh()->Node(me.m_node[k]).get(m_dofP);
                        for (int isol=0; isol<nsol; ++isol) {
                            cm[isol][k] = ms.GetMesh()->Node(me.m_node[k]).get(m_dofC + m_sid[isol]);
                        }
                    }
                    
                    // copy the LM vector
                    ms.UnpackLM(me, mLM);
                    
                    int ndpn;	// number of dofs per node
                    int ndof;	// number of dofs in stiffness matrix
                    
                    if (nsol) {
                        // calculate dofs for biphasic-solute contact
                        ndpn = 4+nsol;
                        ndof = ndpn*(nseln+nmeln);
                        
                        // build the LM vector
                        LM.resize(ndof);
                        
                        for (k=0; k<nseln; ++k)
                        {
                            LM[ndpn*k  ] = sLM[3*k  ];			// x-dof
                            LM[ndpn*k+1] = sLM[3*k+1];			// y-dof
                            LM[ndpn*k+2] = sLM[3*k+2];			// z-dof
                            LM[ndpn*k+3] = sLM[3*nseln+k];		// p-dof
                            for (int isol=0; isol<nsol; ++isol)
                                LM[ndpn*k+4+isol] = sLM[(4+m_sid[isol])*nseln+k];		// c-dof
                        }
                        for (k=0; k<nmeln; ++k)
                        {
                            LM[ndpn*(k+nseln)  ] = mLM[3*k  ];			// x-dof
                            LM[ndpn*(k+nseln)+1] = mLM[3*k+1];			// y-dof
                            LM[ndpn*(k+nseln)+2] = mLM[3*k+2];			// z-dof
                            LM[ndpn*(k+nseln)+3] = mLM[3*nmeln+k];		// p-dof
                            for (int isol=0; isol<nsol; ++isol)
                                LM[ndpn*(k+nseln)+4+isol] = mLM[(4+m_sid[isol])*nmeln+k];		// c-dof
                        }
                    }
                    
                    else if (sporo && mporo) {
                        // calculate dofs for biphasic contact
                        ndpn = 4;
                        ndof = ndpn*(nseln+nmeln);
                        
                        // build the LM vector
                        LM.resize(ndof);
                        
                        for (k=0; k<nseln; ++k)
                        {
                            LM[ndpn*k  ] = sLM[3*k  ];			// x-dof
                            LM[ndpn*k+1] = sLM[3*k+1];			// y-dof
                            LM[ndpn*k+2] = sLM[3*k+2];			// z-dof
                            LM[ndpn*k+3] = sLM[3*nseln+k];		// p-dof
                        }
                        for (k=0; k<nmeln; ++k)
                        {
                            LM[ndpn*(k+nseln)  ] = mLM[3*k  ];			// x-dof
                            LM[ndpn*(k+nseln)+1] = mLM[3*k+1];			// y-dof
                            LM[ndpn*(k+nseln)+2] = mLM[3*k+2];			// z-dof
                            LM[ndpn*(k+nseln)+3] = mLM[3*nmeln+k];		// p-dof
                        }
                    }
                    
                    else {
                        // calculate dofs for elastic contact
                        ndpn = 3;
                        ndof = ndpn*(nseln + nmeln);
                        
                        // build the LM vector
                        LM.resize(ndof);
                        
                        for (k=0; k<nseln; ++k)
                        {
                            LM[3*k  ] = sLM[3*k  ];
                            LM[3*k+1] = sLM[3*k+1];
                            LM[3*k+2] = sLM[3*k+2];
                        }
                        
                        for (k=0; k<nmeln; ++k)
                        {
                            LM[3*(k+nseln)  ] = mLM[3*k  ];
                            LM[3*(k+nseln)+1] = mLM[3*k+1];
                            LM[3*(k+nseln)+2] = mLM[3*k+2];
                        }
                    }
                    
                    // build the en vector
                    en.resize(nseln+nmeln);
                    for (k=0; k<nseln; ++k) en[k      ] = se.m_node[k];
                    for (k=0; k<nmeln; ++k) en[k+nseln] = me.m_node[k];
                    
                    // shape functions
                    Hs = se.H(j);
                    
                    // secondary surface element shape functions
                    double r = pt.m_rs[0];
                    double s = pt.m_rs[1];
                    me.shape_fnc(Hm, r, s);
                    
                    // get normal vector
                    vec3d nu = pt.m_nu;
                    
                    // penalty
                    double eps = m_epsn*pt.m_epsn;
                    
                    // create the stiffness matrix
                    ke.resize(ndof, ndof); ke.zero();
                    
                    // --- S O L I D - S O L I D   C O N T A C T ---
                    
                    mat3ds Ns = mat3dd(1);
                    
                    double* Gr = se.Gr(j);
                    double* Gs = se.Gs(j);
                    vec3d gs[2];
                    ss.CoBaseVectors(se, j, gs);
                    vec3d as[MN];
                    mat3d As[MN];
                    for (int c=0; c<nseln; ++c) {
                        as[c] = (nu ^ (gs[1]*Gr[c] - gs[0]*Gs[c]))/detJ[j];
                        As[c] = (tn[j] & as[c]);
                    }
                    
                    for (int a=0; a<nseln; ++a) {
                        k = a*ndpn;
                        for (int c=0; c<nseln; ++c) {
                            l = c*ndpn;
                            mat3d Kac = (Ns*(Hs[c]*eps)+As[c])*(Hs[a]*detJ[j]*w[j]);
                            ke[k  ][l  ] += Kac(0,0); ke[k  ][l+1] += Kac(0,1); ke[k  ][l+2] += Kac(0,2);
                            ke[k+1][l  ] += Kac(1,0); ke[k+1][l+1] += Kac(1,1); ke[k+1][l+2] += Kac(1,2);
                            ke[k+2][l  ] += Kac(2,0); ke[k+2][l+1] += Kac(2,1); ke[k+2][l+2] += Kac(2,2);
                        }
                        for (int d=0; d<nmeln; ++d) {
                            l = (nseln+d)*ndpn;
                            mat3d Kad = (Ns*(-Hs[a]*Hm[d]*eps))*(detJ[j]*w[j]);
                            ke[k  ][l  ] += Kad(0,0); ke[k  ][l+1] += Kad(0,1); ke[k  ][l+2] += Kad(0,2);
                            ke[k+1][l  ] += Kad(1,0); ke[k+1][l+1] += Kad(1,1); ke[k+1][l+2] += Kad(1,2);
                            ke[k+2][l  ] += Kad(2,0); ke[k+2][l+1] += Kad(2,1); ke[k+2][l+2] += Kad(2,2);
                        }
                    }
                    for (int b=0; b<nmeln; ++b) {
                        k = (nseln+b)*ndpn;
                        for (int c=0; c<nseln; ++c) {
                            l = c*ndpn;
                            mat3d Kbc = (Ns*(Hs[c]*eps)+As[c])*(-Hm[b]*detJ[j]*w[j]);
                            ke[k  ][l  ] += Kbc(0,0); ke[k  ][l+1] += Kbc(0,1); ke[k  ][l+2] += Kbc(0,2);
                            ke[k+1][l  ] += Kbc(1,0); ke[k+1][l+1] += Kbc(1,1); ke[k+1][l+2] += Kbc(1,2);
                            ke[k+2][l  ] += Kbc(2,0); ke[k+2][l+1] += Kbc(2,1); ke[k+2][l+2] += Kbc(2,2);
                        }
                        for (int d=0; d<nmeln; ++d) {
                            l = (nseln+d)*ndpn;
                            mat3d Kbd = Ns*(Hm[b]*Hm[d]*eps*detJ[j]*w[j]);
                            ke[k  ][l  ] += Kbd(0,0); ke[k  ][l+1] += Kbd(0,1); ke[k  ][l+2] += Kbd(0,2);
                            ke[k+1][l  ] += Kbd(1,0); ke[k+1][l+1] += Kbd(1,1); ke[k+1][l+2] += Kbd(1,2);
                            ke[k+2][l  ] += Kbd(2,0); ke[k+2][l+1] += Kbd(2,1); ke[k+2][l+2] += Kbd(2,2);
                        }
                    }
                    
                    // --- M U L T I P H A S I C   S T I F F N E S S ---
                    if (sporo && mporo)
                    {
                        double dt = fem.GetTime().timeIncrement;
                        
                        double epsp = m_epsp*pt.m_epsp;
                        
                        // --- S O L I D - P R E S S U R E / S O L U T E   C O N T A C T ---
                        
                        for (int a=0; a<nseln; ++a) {
                            k = a*ndpn;
                            for (int c=0; c<nseln; ++c) {
                                l = c*ndpn;
                                vec3d gac = (as[c]*(Hs[a]*wn[j]))*(dt*detJ[j]*w[j]);
                                ke[k+3][l  ] += gac.x; ke[k+3][l+1] += gac.y; ke[k+3][l+2] += gac.z;
                                for (int isol=0; isol<nsol; ++isol) {
                                    double epsc = m_epsc*pt.m_epsc[sl[isol]];
                                    vec3d hac = (as[c]*(Hs[a]*jn[isol][j])*epsc)*(dt*detJ[j]*w[j]);
                                    ke[k+4+isol][l  ] += hac.x; ke[k+4+isol][l+1] += hac.y; ke[k+4+isol][l+2] += hac.z;
                                }
                            }
                        }
                        for (int b=0; b<nmeln; ++b) {
                            k = (nseln+b)*ndpn;
                            for (int c=0; c<nseln; ++c) {
                                l = c*ndpn;
                                vec3d gbc = (-as[c]*(Hm[b]*wn[j]))*(dt*detJ[j]*w[j]);
                                ke[k+3][l  ] += gbc.x; ke[k+3][l+1] += gbc.y; ke[k+3][l+2] += gbc.z;
                                for (int isol=0; isol<nsol; ++isol) {
                                    double epsc = m_epsc*pt.m_epsc[sl[isol]];
                                    vec3d hbc = (-as[c]*(Hm[b]*jn[isol][j])*epsc)*(dt*detJ[j]*w[j]);
                                    ke[k+4+isol][l  ] += hbc.x; ke[k+4+isol][l+1] += hbc.y; ke[k+4+isol][l+2] += hbc.z;
                                }
                            }
                        }
                        
                        // --- P R E S S U R E - P R E S S U R E  /  C O N C E N T R A T I O N - C O N C E N T R A T I O N  C O N T A C T ---
                        
                        for (int a=0; a<nseln; ++a) {
                            k = a*ndpn;
                            for (int c=0; c<nseln; ++c) {
                                l = c*ndpn;
                                double gac = (-Hs[a]*Hs[c]*epsp)*(dt*detJ[j]*w[j]);
                                ke[k+3][l+3] += gac;
                                for (int isol=0; isol<nsol; ++isol) {
                                    double epsc = m_epsc*pt.m_epsc[sl[isol]];
                                    double hac = (-Hs[a]*Hs[c]*epsc)*(dt*detJ[j]*w[j]);
                                    ke[k+4+isol][l+4+isol] += hac;
                                }
                            }
                            for (int d=0; d<nmeln; ++d) {
                                l = (nseln+d)*ndpn;
                                double gad = (Hs[a]*Hm[d]*epsp)*(dt*detJ[j]*w[j]);
                                ke[k+3][l+3] += gad;
                                for (int isol=0; isol<nsol; ++isol) {
                                    double epsc = m_epsc*pt.m_epsc[sl[isol]];
                                    double had = (Hs[a]*Hm[d]*epsc)*(dt*detJ[j]*w[j]);
                                    ke[k+4+isol][l+4+isol] += had;
                                }
                            }
                        }
                        for (int b=0; b<nmeln; ++b) {
                            k = (nseln+b)*ndpn;
                            for (int c=0; c<nseln; ++c) {
                                l = c*ndpn;
                                double gbc = (Hm[b]*Hs[c]*epsp)*(dt*detJ[j]*w[j]);
                                ke[k+3][l+3] += gbc;
                                for (int isol=0; isol<nsol; ++isol) {
                                    double epsc = m_epsc*pt.m_epsc[sl[isol]];
                                    double hbc = (Hm[b]*Hs[c]*epsc)*(dt*detJ[j]*w[j]);
                                    ke[k+4+isol][l+4+isol] += hbc;
                                }
                            }
                            for (int d=0; d<nmeln; ++d) {
                                l = (nseln+d)*ndpn;
                                double gbd = (-Hm[b]*Hm[d]*epsp)*(dt*detJ[j]*w[j]);
                                ke[k+3][l+3] += gbd;
                                for (int isol=0; isol<nsol; ++isol) {
                                    double epsc = m_epsc*pt.m_epsc[sl[isol]];
                                    double hbd = (-Hm[b]*Hm[d]*epsc)*(dt*detJ[j]*w[j]);
                                    ke[k+4+isol][l+4+isol] += hbd;
                                }
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
bool FETiedMultiphasicInterface::Augment(int naug, const FETimeInfo& tp)
{
    // make sure we need to augment
	if (m_laugon != 1) return true;

    vec3d Ln;
    double Lp;
    int nsol = (int)m_sid.size();
    vector<double>Lc(nsol);
    bool bconv = true;
    
    bool bporo = (m_ss.m_bporo && m_ms.m_bporo);
    bool bsolu = (m_ss.m_bsolu && m_ms.m_bsolu);

	int NS = m_ss.Elements();
	int NM = m_ms.Elements();

    // --- c a l c u l a t e   i n i t i a l   n o r m s ---
    // a. normal component
    double normL0 = 0, normP = 0, normDP = 0, normC = 0;
    vector<double>normDC(nsol,0);
    for (int i=0; i<NS; ++i)
    {
		FESurfaceElement& el = m_ss.Element(i);
		for (int j = 0; j<el.GaussPoints(); ++j)
		{
			FETiedMultiphasicSurface::Data& ds = static_cast<FETiedMultiphasicSurface::Data&>(*el.GetMaterialPoint(j));
			normL0 += ds.m_Lmd*ds.m_Lmd;
        }
    }
    for (int i=0; i<NM; ++i)
    {
		FESurfaceElement& el = m_ms.Element(i);
        for (int j=0; j<el.GaussPoints(); ++j)
        {
			FETiedMultiphasicSurface::Data& dm = static_cast<FETiedMultiphasicSurface::Data&>(*el.GetMaterialPoint(j));
			normL0 += dm.m_Lmd*dm.m_Lmd;
        }
    }
    
    // b. gap component
    // (is calculated during update)
    double maxgap = 0;
    double maxpg = 0;
    vector<double> maxcg(nsol,0);
    
    // update Lagrange multipliers
    double normL1 = 0, eps, epsp, epsc;
	for (int i = 0; i<NS; ++i)
	{
		FESurfaceElement& el = m_ss.Element(i);
		for (int j = 0; j<el.GaussPoints(); ++j)
		{
			FETiedMultiphasicSurface::Data& ds = static_cast<FETiedMultiphasicSurface::Data&>(*el.GetMaterialPoint(j));

            // update Lagrange multipliers on primary surface
            eps = m_epsn*ds.m_epsn;
            ds.m_Lmd = ds.m_Lmd + ds.m_dg*eps;
            
            normL1 += ds.m_Lmd*ds.m_Lmd;
            
            if (m_ss.m_bporo) {
                Lp = 0;
                Lc.assign(nsol, 0);
                if (ds.m_pme) {
                    epsp = m_epsp*ds.m_epsp;
                    Lp = ds.m_Lmp + epsp*ds.m_pg;
                    maxpg = max(maxpg,fabs(ds.m_pg));
                    normDP += ds.m_pg*ds.m_pg;
                    for (int isol=0; isol<nsol; ++isol) {
                        int l = m_ssl[isol];
                        epsc = m_epsc*ds.m_epsc[l];
                        Lc[isol] = ds.m_Lmc[l] + epsc*ds.m_cg[l];
                        maxcg[isol] = max(maxcg[isol],fabs(ds.m_cg[l]));
                        normDC[isol] += ds.m_cg[l]*ds.m_cg[l];
                    }
                }
                ds.m_Lmp = Lp;
                for (int isol=0; isol<nsol; ++isol)
                    ds.m_Lmc[m_ssl[isol]] = Lc[isol];
            }
            
            maxgap = max(maxgap,sqrt(ds.m_dg*ds.m_dg));
        }
    }
    
	for (int i = 0; i<NM; ++i)
	{
		FESurfaceElement& el = m_ms.Element(i);
		for (int j = 0; j<el.GaussPoints(); ++j)
		{
			FETiedMultiphasicSurface::Data& dm = static_cast<FETiedMultiphasicSurface::Data&>(*el.GetMaterialPoint(j));

            // update Lagrange multipliers on secondary surface
            eps = m_epsn*dm.m_epsn;
            dm.m_Lmd = dm.m_Lmd + dm.m_dg*eps;
            
            normL1 += dm.m_Lmd*dm.m_Lmd;
            
            if (m_ms.m_bporo) {
                Lp = 0;
                Lc.assign(nsol, 0);
                if (dm.m_pme) {
                    epsp = m_epsp*dm.m_epsp;
                    Lp = dm.m_Lmp + epsp*dm.m_pg;
                    maxpg = max(maxpg,fabs(dm.m_pg));
                    normDP += dm.m_pg*dm.m_pg;
                    for (int isol=0; isol<nsol; ++isol) {
                        int l = m_msl[isol];
                        epsc = m_epsc*dm.m_epsc[l];
                        Lc[isol] = dm.m_Lmc[l] + epsc*dm.m_cg[l];
                        maxcg[isol] = max(maxcg[isol],fabs(dm.m_cg[l]));
                        normDC[isol] += dm.m_cg[l]*dm.m_cg[l];
                    }
                }
                dm.m_Lmp = Lp;
                for (int isol=0; isol<nsol; ++isol)
                    dm.m_Lmc[m_msl[isol]] = Lc[isol];
            }
            
            maxgap = max(maxgap,sqrt(dm.m_dg*dm.m_dg));
        }
    }
    
    // Ideally normP should be evaluated from the fluid pressure at the
    // contact interface (not easily accessible).  The next best thing
    // is to use the contact traction.
    normP = normL1;
    normC = normL1/(m_Rgas*m_Tabs);
    
    // calculate relative norms
    double lnorm = (normL1 != 0 ? fabs((normL1 - normL0) / normL1) : fabs(normL1 - normL0)); 
    double pnorm = (normP != 0 ? (normDP/normP) : normDP); 
    vector<double> cnorm(nsol);
    for (int isol=0; isol<nsol; ++isol)
        cnorm[isol] = (normC != 0 ? (normDC[isol]/normC) : normDC[isol]);
    
    // check convergence
    if ((m_gtol > 0) && (maxgap > m_gtol)) bconv = false;
    if ((m_ptol > 0) && (bporo && maxpg > m_ptol)) bconv = false;
    for (int isol=0; isol<nsol; ++isol)
        if ((m_ctol > 0) && (bsolu && maxcg[isol] > m_ctol)) bconv = false;
    
    if ((m_atol > 0) && (lnorm > m_atol)) bconv = false;
    if ((m_atol > 0) && (pnorm > m_atol)) bconv = false;
    for (int isol=0; isol<nsol; ++isol)
        if ((m_atol > 0) && (cnorm[isol] > m_atol)) bconv = false;
    
    if (naug < m_naugmin ) bconv = false;
    if (naug >= m_naugmax) bconv = true;
    
    feLog(" sliding interface # %d\n", GetID());
    feLog("                        CURRENT        REQUIRED\n");
    feLog("    D multiplier : %15le", lnorm);
    if (m_atol > 0) feLog("%15le\n", m_atol);
    else feLog("       ***\n");
    if (bporo) {
        feLog("    P gap        : %15le", pnorm);
        if (m_atol > 0) feLog("%15le\n", m_atol);
        else feLog("       ***\n"); }
    for (int isol=0; isol<nsol; ++isol) {
        feLog("    C[%d] gap   : %15le", m_sid[isol], cnorm[isol]);
        if (m_atol > 0) feLog("%15le\n", m_atol);
        else feLog("       ***\n");
    }
    
    feLog("    maximum gap  : %15le", maxgap);
    if (m_gtol > 0) feLog("%15le\n", m_gtol);
    else feLog("       ***\n");
    if (bporo) {
        feLog("    maximum pgap : %15le", maxpg);
        if (m_ptol > 0) feLog("%15le\n", m_ptol);
        else feLog("       ***\n");
    }
    for (int isol=0; isol<nsol; ++isol) {
        feLog("    maximum cgap[%d] : %15le", m_sid[isol], maxcg[isol]);
        if (m_ctol > 0) feLog("%15le\n", m_ctol);
        else feLog("       ***\n");
    }
    
    return bconv;
}

//-----------------------------------------------------------------------------
void FETiedMultiphasicInterface::Serialize(DumpStream &ar)
{
    // store base class data
    FEContactInterface::Serialize(ar);
    
    // store contact surface data
    m_ms.Serialize(ar);
    m_ss.Serialize(ar);
    
	// serialize interface data
	ar & m_epsp;
	ar & m_epsc;
	ar & m_sid;
	ar & m_ssl;
	ar & m_msl;
	ar & m_sz;

	// serialize element pointers
	SerializeElementPointers(m_ss, m_ms, ar);
	SerializeElementPointers(m_ms, m_ss, ar);
}

//-----------------------------------------------------------------------------
FESoluteData* FETiedMultiphasicInterface::FindSoluteData(int nid)
{
    FEModel& fem = *GetFEModel();
    int N = GetFEModel()->GlobalDataItems();
    for (int i=0; i<N; ++i)
    {
        FESoluteData* psd = dynamic_cast<FESoluteData*>(fem.GetGlobalData(i));
        if (psd && (psd->GetID() == nid)) return psd;
    }
    return 0;
}
