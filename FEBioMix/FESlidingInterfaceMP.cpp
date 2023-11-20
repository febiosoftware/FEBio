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
#include "FESlidingInterfaceMP.h"
#include "FEBiphasic.h"
#include "FEBiphasicSolute.h"
#include "FETriphasic.h"
#include "FEMultiphasic.h"
#include "FECore/FEModel.h"
#include "FECore/log.h"
#include "FECore/DOFS.h"
#include "FECore/FENormalProjection.h"
#include "FECore/FEAnalysis.h"
#include <FECore/FELinearSystem.h>

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEAmbientConcentration, FECoreClass)
	ADD_PARAMETER(m_sol, "sol", FE_PARAM_ATTRIBUTE, "$(solutes)");
	ADD_PARAMETER(m_ambc, "ambient_concentration")->MakeVolatile(false);
END_FECORE_CLASS();

FEAmbientConcentration::FEAmbientConcentration(FEModel* fem) : FECoreClass(fem) 
{
	m_sol = -1;
	m_ambc = 0.0;
}

//-----------------------------------------------------------------------------
// Define sliding interface parameters
BEGIN_FECORE_CLASS(FESlidingInterfaceMP, FEContactInterface)
    ADD_PARAMETER(m_atol     , "tolerance"            );
    ADD_PARAMETER(m_gtol     , "gaptol"               )->setUnits(UNIT_LENGTH);;
    ADD_PARAMETER(m_ptol     , "ptol"                 );
    ADD_PARAMETER(m_ctol     , "ctol"                 );
    ADD_PARAMETER(m_epsn     , "penalty"              );
    ADD_PARAMETER(m_bautopen , "auto_penalty"         );
    ADD_PARAMETER(m_bupdtpen , "update_penalty"       );
    ADD_PARAMETER(m_btwo_pass, "two_pass"             );
    ADD_PARAMETER(m_knmult   , "knmult"               );
    ADD_PARAMETER(m_stol     , "search_tol"           );
    ADD_PARAMETER(m_epsp     , "pressure_penalty"     );
    ADD_PARAMETER(m_epsc     , "concentration_penalty");
    ADD_PARAMETER(m_bsymm    , "symmetric_stiffness"  );
    ADD_PARAMETER(m_srad     , "search_radius"        )->setUnits(UNIT_LENGTH);;
    ADD_PARAMETER(m_nsegup   , "seg_up"               );
    ADD_PARAMETER(m_breloc   , "node_reloc"           );
    ADD_PARAMETER(m_mu       , "fric_coeff"           );
    ADD_PARAMETER(m_phi      , "contact_frac"         );
    ADD_PARAMETER(m_bsmaug   , "smooth_aug"           );
    ADD_PARAMETER(m_bsmfls   , "smooth_fls"           );
    ADD_PARAMETER(m_naugmin  , "minaug"               );
    ADD_PARAMETER(m_naugmax  , "maxaug"               );
    ADD_PARAMETER(m_ambp     , "ambient_pressure"     );

	ADD_PROPERTY(m_ambctmp, "ambient_concentration",FEProperty::Optional);
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FESlidingSurfaceMP::Data::Data()
{
    m_Lmd   = 0.0;
    m_Lmt   = m_tr = vec3d(0,0,0);
    m_Lmp   = 0.0;
    m_epsn  = 1.0;
    m_epsp  = 1.0;
    m_pg    = 0.0;
    m_p1    = 0.0;
    m_mueff = 0.0;
    m_fls   = 0.0;
    m_nu    = m_s1 = m_dg = vec3d(0,0,0);
    m_rs    = m_rsp = vec2d(0,0);
    m_bstick = false;
}

//-----------------------------------------------------------------------------
void FESlidingSurfaceMP::Data::Init()
{
    m_Lmd   = 0.0;
    m_Lmt   = m_tr = vec3d(0,0,0);
    m_Lmp   = 0.0;
    m_epsn  = 1.0;
    m_epsp  = 1.0;
    m_pg    = 0.0;
    m_p1    = 0.0;
    m_mueff = 0.0;
    m_fls   = 0.0;
    m_nu    = m_s1 = m_dg = vec3d(0,0,0);
    m_rs    = m_rsp = vec2d(0,0);
    m_bstick = false;
}

//-----------------------------------------------------------------------------
void FESlidingSurfaceMP::Data::Serialize(DumpStream& ar)
{
	FEBiphasicContactPoint::Serialize(ar);
    ar & m_gap;
    ar & m_dg;
    ar & m_nu;
    ar & m_s1;
    ar & m_rs;
    ar & m_rsp;
    ar & m_Lmd;
    ar & m_Lmt;
    ar & m_Lmp;
    ar & m_Lmc;
    ar & m_epsn;
    ar & m_epsp;
    ar & m_epsc;
    ar & m_pg;
    ar & m_cg;
    ar & m_p1;
    ar & m_c1;
    ar & m_Ln;
    ar & m_bstick;
    ar & m_tr;
    ar & m_mueff;
    ar & m_fls;
}

//-----------------------------------------------------------------------------
// FESlidingSurfaceMP
//-----------------------------------------------------------------------------

FESlidingSurfaceMP::FESlidingSurfaceMP(FEModel* pfem) : FEBiphasicContactSurface(pfem)
{ 
	m_bporo = m_bsolu = false;
	m_dofC = -1;
}

//-----------------------------------------------------------------------------
FESlidingSurfaceMP::~FESlidingSurfaceMP()
{

}

//-----------------------------------------------------------------------------
//! create material point data
FEMaterialPoint* FESlidingSurfaceMP::CreateMaterialPoint()
{
	return new FESlidingSurfaceMP::Data;
}

//-----------------------------------------------------------------------------
void FESlidingSurfaceMP::UnpackLM(FEElement& el, vector<int>& lm)
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
bool FESlidingSurfaceMP::Init()
{
	// initialize surface data first
	if (FEBiphasicContactSurface::Init() == false) return false;

	// store concentration index
	DOFS& dofs = GetFEModel()->GetDOFS();
	m_dofC = dofs.GetDOF("concentration", 0);
	
    // allocate node normals and pressures
    m_nn.assign(Nodes(), vec3d(0,0,0));
    m_tn.assign(Nodes(), vec3d(0,0,0));
    m_pn.assign(Nodes(), 0);

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
            FETriphasic* ptp = dynamic_cast<FETriphasic*> (pm);
			FEMultiphasic* pmp = dynamic_cast<FEMultiphasic*> (pm);
			if (pb) {
				m_bporo = true;
				nsol = 0;
			}
            else if (pbs) {
				m_bporo = m_bsolu = true;
				nsol = 1;
				m_sid.assign(nsol, pbs->GetSolute()->GetSoluteID() - 1);
			}
            else if (ptp) {
                m_bporo = m_bsolu = true;
                nsol = ptp->Solutes();
                m_sid.resize(nsol);
                for (int isol=0; isol<nsol; ++isol) {
                    m_sid[isol] = ptp->GetSolute(isol)->GetSoluteID() - 1;
                }
            }
            else if (pmp) {
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
	for (int i=0; i<NE; ++i)
	{
		FESurfaceElement& el = Element(i);
		int nint = el.GaussPoints();
        if (nsol) {
            for (int j=0; j<nint; ++j) {
				Data& data = static_cast<Data&>(*el.GetMaterialPoint(j));
                data.m_Lmc.resize(nsol);
                data.m_epsc.resize(nsol);
                data.m_cg.resize(nsol);
                data.m_c1.resize(nsol);
                data.m_epsc.assign(nsol,1);
                data.m_c1.assign(nsol,0.0);
                data.m_cg.assign(nsol,0.0);
                data.m_Lmc.assign(nsol,0.0);
            }
        }
	}
	
	return true;
}

//-----------------------------------------------------------------------------
void FESlidingSurfaceMP::InitSlidingSurface()
{
    for (int i=0; i<Elements(); ++i)
    {
        FESurfaceElement& el = Element(i);
        int nint = el.GaussPoints();
        for (int j=0; j<nint; ++j)
        {
            // Store current surface projection values as previous
            Data& data = static_cast<Data&>(*el.GetMaterialPoint(j));
            data.m_rsp = data.m_rs;
            data.m_pmep = data.m_pme;
        }
    }
}

//-----------------------------------------------------------------------------
//! Evaluate the nodal contact pressures by averaging values from surrounding
//! faces.  This function ensures that nodal contact pressures are always
//! positive, so that they can be used to detect free-draining status.

void FESlidingSurfaceMP::EvaluateNodalContactPressures()
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
//! Evaluate the nodal contact pressures by averaging values from surrounding
//! faces.  This function ensures that nodal contact pressures are always
//! positive, so that they can be used to detect free-draining status.

void FESlidingSurfaceMP::EvaluateNodalContactTractions()
{
    const int N = Nodes();
    
    // number of faces with non-zero contact pressure connected to this node
    vector<int> nfaces(N,0);
    
    // zero nodal contact pressures
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

void FESlidingSurfaceMP::UpdateNodeNormals()
{
	int N = Nodes(), i, j, ne, jp1, jm1;
	const int MN = FEElement::MAX_NODES;
	vec3d y[MN], n;
	
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
vec3d FESlidingSurfaceMP::GetContactForce()
{
    return m_Ft;
}

//-----------------------------------------------------------------------------
double FESlidingSurfaceMP::GetContactArea()
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
			double s = (data.m_Ln > 0) ? 1 : 0;
            
			// get the base vectors
			vec3d g[2];
			CoBaseVectors(el, i, g);
            
			// normal (magnitude = area)
			vec3d n = g[0] ^ g[1];
            
			// gauss weight
			double w = el.GaussWeights()[i];
            
			// contact force
			a += n.norm()*(w*s);
		}
	}
	
	return a;
}

//-----------------------------------------------------------------------------
vec3d FESlidingSurfaceMP::GetFluidForce()
{
	int n, i;
	const int MN = FEElement::MAX_NODES;
	double pn[MN];
	
    // get parent contact interface to extract ambient fluid pressure
    FESlidingInterfaceMP* simp = dynamic_cast<FESlidingInterfaceMP*>(GetContactInterface());
    double ambp = simp->m_ambp;
    
	// initialize contact force
	vec3d f(0,0,0);
	
	// loop over all elements of the surface
	for (n=0; n<Elements(); ++n)
	{
		FESurfaceElement& el = Element(n);
		int nseln = el.Nodes();
		
		// nodal pressures
		for (i=0; i<nseln; ++i) pn[i] = GetMesh()->Node(el.m_node[i]).get(m_dofP);
		
		int nint = el.GaussPoints();
		
		// evaluate the fluid force for that element
		for (i=0; i<nint; ++i) 
		{
			// get the base vectors
			vec3d g[2];
			CoBaseVectors(el, i, g);
			// normal (magnitude = area)
			vec3d n = g[0] ^ g[1];
			// gauss weight
			double w = el.GaussWeights()[i];
            // fluid pressure for fluid load support
            double p = el.eval(pn, i) - ambp;
			// contact force
			f += n*(w*p);
		}
	}
	
	return f;
}

//-----------------------------------------------------------------------------
void FESlidingSurfaceMP::Serialize(DumpStream& ar)
{
	FEBiphasicContactSurface::Serialize(ar);
	ar & m_dofP & m_dofC;
	ar & m_bporo;
	ar & m_bsolu;
	ar & m_nn;
	ar & m_sid;
	ar & m_pn;
    ar & m_tn;
	ar & m_Ft;
}

//-----------------------------------------------------------------------------
void FESlidingSurfaceMP::GetVectorGap(int nface, vec3d& pg)
{
    FESurfaceElement& el = Element(nface);
    int ni = el.GaussPoints();
    pg = vec3d(0,0,0);
    for (int k=0; k<ni; ++k)
    {
        Data& data = static_cast<Data&>(*el.GetMaterialPoint(k));
        pg += data.m_dg;
    }
    pg /= ni;
}

//-----------------------------------------------------------------------------
void FESlidingSurfaceMP::GetContactPressure(int nface, double& pg)
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
void FESlidingSurfaceMP::GetContactTraction(int nface, vec3d& pt)
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
void FESlidingSurfaceMP::GetSlipTangent(int nface, vec3d& pt)
{
    FESurfaceElement& el = Element(nface);
    int ni = el.GaussPoints();
    pt = vec3d(0,0,0);
    for (int k=0; k<ni; ++k)
    {
        Data& data = static_cast<Data&>(*el.GetMaterialPoint(k));
        if (!data.m_bstick) pt += data.m_s1;
    }
    pt /= ni;
}

//-----------------------------------------------------------------------------
void FESlidingSurfaceMP::GetMuEffective(int nface, double& pg)
{
    FESurfaceElement& el = Element(nface);
    int ni = el.GaussPoints();
    pg = 0;
    for (int k=0; k<ni; ++k)
    {
        Data& data = static_cast<Data&>(*el.GetMaterialPoint(k));
        pg += data.m_mueff;
    }
    pg /= ni;
}

//-----------------------------------------------------------------------------
void FESlidingSurfaceMP::GetLocalFLS(int nface, double& pg)
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
void FESlidingSurfaceMP::GetNodalVectorGap(int nface, vec3d* pg)
{
    FESurfaceElement& el = Element(nface);
    int ni = el.GaussPoints();
    vec3d gi[FEElement::MAX_INTPOINTS];
    for (int k=0; k<ni; ++k)
    {
        Data& data = static_cast<Data&>(*el.GetMaterialPoint(k));
        gi[k] = data.m_dg;
    }
    el.project_to_nodes(gi, pg);
}

//-----------------------------------------------------------------------------
void FESlidingSurfaceMP::GetNodalContactPressure(int nface, double* pn)
{
	FESurfaceElement& el = Element(nface);
    for (int k=0; k<el.Nodes(); ++k)
        pn[k] = m_pn[el.m_lnode[k]];
}

//-----------------------------------------------------------------------------
void FESlidingSurfaceMP::GetStickStatus(int nface, double& pg)
{
    FESurfaceElement& el = Element(nface);
    int ni = el.GaussPoints();
    pg = 0;
    for (int k=0; k<ni; ++k)
    {
        Data& data = static_cast<Data&>(*el.GetMaterialPoint(k));
        if (data.m_bstick) pg += 1.0;
    }
    pg /= ni;
}

//-----------------------------------------------------------------------------
void FESlidingSurfaceMP::GetNodalContactTraction(int nface, vec3d* tn)
{
	FESurfaceElement& el = Element(nface);
    for (int k=0; k<el.Nodes(); ++k)
        tn[k] = m_tn[el.m_lnode[k]];
}

//-----------------------------------------------------------------------------
// FESlidingInterfaceMP
//-----------------------------------------------------------------------------

FESlidingInterfaceMP::FESlidingInterfaceMP(FEModel* pfem) : FEContactInterface(pfem), m_ss(pfem), m_ms(pfem)
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
	m_ambp = 0;
	m_nsegup = 0;
	m_bautopen = false;
    m_breloc = false;
    m_bsmaug = false;
    m_bupdtpen = false;
    m_mu = 0.0;
    m_phi = 0.0;

	m_naugmin = 0;
	m_naugmax = 10;

    m_bfreeze = false;

	m_dofP = -1;
	m_dofC = -1;

	m_ss.SetSibling(&m_ms);
	m_ms.SetSibling(&m_ss);
    m_ss.SetContactInterface(this);
    m_ms.SetContactInterface(this);
}

//-----------------------------------------------------------------------------

FESlidingInterfaceMP::~FESlidingInterfaceMP()
{
}

//-----------------------------------------------------------------------------
bool FESlidingInterfaceMP::Init()
{
	m_Rgas = GetFEModel()->GetGlobalConstant("R");
	m_Tabs = GetFEModel()->GetGlobalConstant("T");

	// get number of DOFS
	FEModel* fem = GetFEModel();
	DOFS& fedofs = fem->GetDOFS();
	int nsol = fedofs.GetVariableSize("concentration");
	m_ambc.assign(nsol, 0.0);
	m_dofP = fem->GetDOFIndex("p");
	m_dofC = fem->GetDOFIndex("concentration", 0);
	
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
                FESoluteData* sd = FindSoluteData(m_ss.m_sid[is]+1);
                m_sz.push_back(sd->m_z);
			}
		}
	}

    // cycle through all the solutes and determine ambient concentrations
	for (int i = 0; i < m_ambctmp.size(); ++i)
	{
		FEAmbientConcentration* aci = m_ambctmp[i];
		int isol = aci->m_sol - 1;
		assert((isol >= 0) && (isol < nsol));
		m_ambc[isol] = aci->m_ambc;
	}

	return true;
}

//-----------------------------------------------------------------------------
//! build the matrix profile for use in the stiffness matrix
void FESlidingInterfaceMP::BuildMatrixProfile(FEGlobalMatrix& K)
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
		FESlidingSurfaceMP& ss = (np == 0? m_ss : m_ms);
        
		int k, l;
		for (int j=0; j<ss.Elements(); ++j)
		{
			FESurfaceElement& se = ss.Element(j);
			int nint = se.GaussPoints();
			int* sn = &se.m_node[0];
			for (k=0; k<nint; ++k)
			{
				FESlidingSurfaceMP::Data& data = static_cast<FESlidingSurfaceMP::Data&>(*se.GetMaterialPoint(k));

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
void FESlidingInterfaceMP::UpdateAutoPenalty()
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
void FESlidingInterfaceMP::Activate()
{
	// don't forget to call the base members
	FEContactInterface::Activate();
	
    UpdateAutoPenalty();
    
	// update sliding interface data
	Update();
}

//-----------------------------------------------------------------------------
void FESlidingInterfaceMP::CalcAutoPenalty(FESlidingSurfaceMP& s)
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
			FESlidingSurfaceMP::Data& pt = static_cast<FESlidingSurfaceMP::Data&>(*el.GetMaterialPoint(j));
			pt.m_epsn = eps;
        }
	}
}

//-----------------------------------------------------------------------------
void FESlidingInterfaceMP::CalcAutoPressurePenalty(FESlidingSurfaceMP& s)
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
			FESlidingSurfaceMP::Data& pt = static_cast<FESlidingSurfaceMP::Data&>(*el.GetMaterialPoint(j));
			pt.m_epsp = eps;
        }
	}
}

//-----------------------------------------------------------------------------
//! This function calculates a contact penalty parameter based on the
//! material and geometrical properties of the primary and secondary surfaces
//!
double FESlidingInterfaceMP::AutoPenalty(FESurfaceElement& el, FESurface &s)
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
        FETriphasic* pms = dynamic_cast<FETriphasic*>(pme);
        S = pms->Tangent(mp);
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

double FESlidingInterfaceMP::AutoPressurePenalty(FESurfaceElement& el, FESlidingSurfaceMP& s)
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
    if (pb)
        K = pb->GetPermeability()->Permeability(mp);
	else if (pbs)
		K = pbs->GetPermeability()->Permeability(mp);
    else if (ptp)
        K = ptp->GetPermeability()->Permeability(mp);
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
void FESlidingInterfaceMP::CalcAutoConcentrationPenalty(FESlidingSurfaceMP& s,
														const int isol)
{
	// loop over all surface elements
	int ni = 0;
	for (int i=0; i<s.Elements(); ++i)
	{
		// get the surface element
		FESurfaceElement& el = s.Element(i);
		
		// calculate a modulus
		double eps = AutoConcentrationPenalty(el, s, isol);
		
		// assign to integation points of surface element
		int nint = el.GaussPoints();
		for (int j=0; j<nint; ++j, ++ni)
        {
			FESlidingSurfaceMP::Data& pt = static_cast<FESlidingSurfaceMP::Data&>(*el.GetMaterialPoint(j));
			pt.m_epsc[isol] = eps;
        }
	}
}

//-----------------------------------------------------------------------------

double FESlidingInterfaceMP::AutoConcentrationPenalty(FESurfaceElement& el, 
													  FESlidingSurfaceMP& s,
													  const int isol)
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
void FESlidingInterfaceMP::ProjectSurface(FESlidingSurfaceMP& ss, FESlidingSurfaceMP& ms, bool bupseg, bool bmove)
{
    FEMesh& mesh = GetFEModel()->GetMesh();
    FESurfaceElement* pme;
    vec3d r, nu;
    double rs[2] = {0,0};
    double Ln;
    
    const int MN = FEElement::MAX_NODES;
    int nsol = (int)m_sid.size();
    double ps[MN], p1 = 0.0;
    vector< vector<double> > cs(nsol, vector<double>(MN));
    vector<double> c1(nsol);
    c1.assign(nsol,0);
    
    double psf = GetPenaltyScaleFactor();
    
    // initialize projection data
    FENormalProjection np(ms);
    np.SetTolerance(m_stol);
    np.SetSearchRadius(m_srad);
    np.Init();
    
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
    // TODO: commenting this pragma line made my code recover slidingelastic test case results
//#pragma omp parallel for
    for (int i=0; i<ss.Elements(); ++i)
    {
        FESurfaceElement& el = ss.Element(i);
        
        bool sporo = ss.m_bporo;
        
        int ne = el.Nodes();
        int nint = el.GaussPoints();
        
        // get the nodal pressures
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
            FESlidingSurfaceMP::Data& pt = static_cast<FESlidingSurfaceMP::Data&>(*el.GetMaterialPoint(j));

            // calculate the global position of the integration point
            r = ss.Local2Global(el, j);
            
            // get the pressure at the integration point
            if (sporo) p1 = el.eval(ps, j);
            
            // get the concentration at the integration point
            for (int isol=0; isol<nsol; ++isol) c1[isol] = el.eval(&cs[isol][0], j);
            
            // calculate the normal at this integration point
            nu = ss.SurfaceNormal(el, j);
            
            // first see if the old intersected face is still good enough
            pme = pt.m_pme;
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
                
                Ln = pt.m_Lmd + eps*g;
                
                pt.m_gap = (g <= m_srad? g : 0);
                
                bool mporo = ms.m_bporo;
                
                if ((Ln >= 0) && (g <= m_srad))
                {
                    
                    // get the pressure at the contact point
                    // account for mixed multiphasic-elastic contact with elastic primary
                    // calculate the pressure gap function
                    double p2 = 0;
                    if (mporo) {
                        double pm[MN];
                        for (int k=0; k<pme->Nodes(); ++k) pm[k] = mesh.Node(pme->m_node[k]).get(m_dofP);
                        p2 = pme->eval(pm, rs[0], rs[1]);
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

                    for (int isol=0; isol<nsol; ++isol) {
                        int sid = m_sid[isol];
                        double cm[MN];
                        for (int k=0; k<pme->Nodes(); ++k) cm[k] = mesh.Node(pme->m_node[k]).get(m_dofC + sid);
                        double c2 = pme->eval(cm, rs[0], rs[1]);
                        pt.m_cg[m_ssl[isol]] = c1[isol] - c2;
                        pt.m_c1[m_ssl[isol]] = c1[isol];
                    }
                }
                else
                {
                    pt.m_Lmd = 0;
                    pt.m_gap = 0;
                    pt.m_pme = 0;
                    pt.m_dg = pt.m_Lmt = vec3d(0,0,0);
                    if (sporo || mporo) {
                        pt.m_Lmp = 0;
                        pt.m_pg = 0;
                        pt.m_p1 = 0;
                    }
                    for (int isol=0; isol<nsol; ++isol) {
                        pt.m_Lmc[m_ssl[isol]] = 0;
                        pt.m_cg[m_ssl[isol]] = 0;
                        pt.m_c1[m_ssl[isol]] = 0;
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
                for (int isol=0; isol<nsol; ++isol) {
                    pt.m_Lmc[m_ssl[isol]] = 0;
                    pt.m_cg[m_ssl[isol]] = 0;
                    pt.m_c1[m_ssl[isol]] = 0;
                }
            }
        }
    }
}

//-----------------------------------------------------------------------------

void FESlidingInterfaceMP::Update()
{
    double rs[2]={0,0};
    
    FEModel& fem = *GetFEModel();
    
    // get number of DOFS
    DOFS& fedofs = GetFEModel()->GetDOFS();
    int MAX_CDOFS = fedofs.GetVariableSize("concentration");
    
    static int naug = 0;
    static int biter = 0;
    
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
    //    Logfile& log = GetLogfile();
    //    log.printf("seg_up iteration # %d\n", niter+1);
    
    // project the surfaces onto each other
    // this will update the gap functions as well
    static bool bfirst = true;
    ProjectSurface(m_ss, m_ms, bupseg, (m_breloc && bfirst));
    // TODO: there was a bug below - the right part of the OR statement was m_ss.m_bporo
    if (m_btwo_pass || m_ms.m_bporo) ProjectSurface(m_ms, m_ss, bupseg);
    bfirst = false;
    
    // Call InitSlidingSurface on the first iteration of each time step
    // TODO: previously had the line below controlling InitSlidingSurface, but SlidingBiphasic uses nsolve_iter
    //  if (niter == 0)
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

    // TODO: this line is above the "only continue if..." line in SlidingMP
    // update node normals
    m_ss.UpdateNodeNormals();
    m_ms.UpdateNodeNormals();
    
    // Now that the nodes have been projected, we need to figure out
    // if we need to modify the constraints on the pressure and concentration dofs.
    // If the nodes are not in contact, they must match ambient
    // conditions. Since all nodes have been previously marked to be
    // under ambient conditions in MarkAmbient(), we just need to reverse
    // this setting here, for nodes that are in contact.
    
    // Next, we loop over each surface, visiting the nodes
    // and finding out if that node is in contact or not
    int npass = (m_btwo_pass?2:1);
    for (int np=0; np<npass; ++np)
    {
        FESlidingSurfaceMP& ss = (np == 0? m_ss : m_ms);
        FESlidingSurfaceMP& ms = (np == 0? m_ms : m_ss);
        
        // initialize projection data
        FENormalProjection project(ss);
        project.SetTolerance(m_stol);
        project.SetSearchRadius(m_srad);
        project.Init();

        // loop over all the nodes of the primary surface
        for (int n=0; n<ss.Nodes(); ++n) {
            if (ss.m_pn[n] > 0)
            {
                FENode& node = ss.Node(n);
                int id = node.m_ID[m_dofP];
                if (id < -1)
                    // mark node as non-ambient (= pos ID)
                    node.m_ID[m_dofP] = -id-2;
                for (int j=0; j<MAX_CDOFS; ++j) {
                    id = node.m_ID[m_dofC+j];
                    if (id < -1)
                        // mark node as non-ambient (= pos ID)
                        node.m_ID[m_dofC+j] = -id-2;
                }
            }
        }
        
        // loop over all nodes of the secondary surface
        // the secondary surface is trickier since we need
        // to look at the primary's surface projection
        if (ms.m_bporo) {
            for (int n=0; n<ms.Nodes(); ++n)
            {
                // get the node
                FENode& node = ms.Node(n);
                
                // project it onto the primary surface
                FESurfaceElement* pse = project.Project(node.m_rt, ms.m_nn[n], rs);
                
                if (pse)
                {
                    // we found an element, so let's see if it's even remotely close to contact
                    // find the global location of the intersection point
                    vec3d q = ss.Local2Global(*pse, rs[0], rs[1]);
                    
                    // calculate the gap function
                    double g = ms.m_nn[n]*(node.m_rt - q);
                    
                    // TODO: SlidingMP uses R, SlidingBiphasic uses m_srad
                    if (fabs(g) <= m_srad)
                    {
                        // we found an element so let's calculate the nodal traction values for this element
                        // get the normal tractions at the nodes
                        double tn[FEElement::MAX_NODES];
                        for (int i=0; i<pse->Nodes(); ++i)
                            tn[i] = ss.m_pn[pse->m_lnode[i]];
                        
                        // now evaluate the traction at the intersection point
                        double tp = pse->eval(tn, rs[0], rs[1]);
                        
                        // if tp > 0, mark node as non-ambient. (= pos ID)
                        if (tp > 0)  {
                            int id = node.m_ID[m_dofP];
                            if (id < -1)
                                // mark as non-ambient
                                node.m_ID[m_dofP] = -id-2;
                            for (int j=0; j<MAX_CDOFS; ++j) {
                                id = node.m_ID[m_dofC+j];
                                if (id < -1)
                                    // mark as non-ambient
                                    node.m_ID[m_dofC+j] = -id-2;
                            }
                        }
                    }
                }
            }
        }
    }
}

//-----------------------------------------------------------------------------
vec3d FESlidingInterfaceMP::SlipTangent(FESlidingSurfaceMP& ss, const int nel, const int nint, FESlidingSurfaceMP& ms, double& dh, vec3d& r)
{
    vec3d s1(0,0,0);
    dh = 0;
    r = vec3d(0,0,0);
    
    // get primary surface element
    FESurfaceElement& se = ss.Element(nel);
    
    // get integration point data
    FESlidingSurfaceMP::Data& data = static_cast<FESlidingSurfaceMP::Data&>(*se.GetMaterialPoint(nint));
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
vec3d FESlidingInterfaceMP::ContactTraction(FESlidingSurfaceMP& ss, const int nel, const int n, FESlidingSurfaceMP& ms, double& pn)
{
    vec3d s1(0,0,0);
    vec3d dr(0,0,0);
    vec3d t(0,0,0);
    pn = 0;
    double tn = 0, ts = 0, mueff = 0;
    double psf = GetPenaltyScaleFactor();
    
    int nsol = (int)m_sid.size();
    vector<double> c1(nsol);
    c1.assign(nsol,0);
    
    // get the mesh
    FEMesh& m = GetFEModel()->GetMesh();
    
    // get the primary surface element
    FESurfaceElement& se = ss.Element(nel);
    
    // get the integration point data
    FESlidingSurfaceMP::Data& data = static_cast<FESlidingSurfaceMP::Data&>(*se.GetMaterialPoint(n));
    
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
    
    // get the solute concentrations at this integration point
    for (int isol=0; isol<nsol; ++isol) c1[isol] = data.m_c1[m_ssl[isol]];
    
    // get the fluid pressure for load sharing at this integration point
    double ph = data.m_p1 - m_ambp;
    
    // get poro status of primary surface
    bool sporo = ss.m_bporo;
    
    // get current and previous secondary elements
    FESurfaceElement* pme = data.m_pme;
    FESurfaceElement* pmep = data.m_pmep;
    
    // zero the effective friction coefficient
    data.m_mueff = 0.0;
    data.m_fls = 0.0;
    data.m_s1 = vec3d(0,0,0);
    
    // get local FLS from element projection
    // TODO: Gerard: I added an optional 'pamb' argument to GetGPLocalFLS
    double fls = 0;
    if (m_bsmfls) {
        double lfls[FEElement::MAX_INTPOINTS];
        ss.GetGPLocalFLS(nel, lfls, m_ambp);
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
                data.m_fls = m_bsmfls ? fls : ph/pn;
            }
            
            // store the previous values as the current
            data.m_pme = data.m_pmep;
            data.m_rs = data.m_rsp;
            
            // recalculate gap
            data.m_dg = dg;
            
            // recalculate pressure gap
            bool mporo = ms.m_bporo;
            if (sporo && mporo)
            {
                double pm[FEElement::MAX_NODES];
                for (int k=0; k<pme->Nodes(); ++k) pm[k] = m.Node(pme->m_node[k]).get(m_dofP);
                double p2 = pme->eval(pm, data.m_rs[0], data.m_rs[1]);
                data.m_pg = p - p2;
            }
            //TODO: what if stick/slip moves from poro to elastic or vice versa? Should I add
            // an "else {data.m_pg = 0}" line?
            
            // recalculate concentration gaps
            for (int isol=0; isol<nsol; ++isol) {
                int sid = m_sid[isol];
                double cm[FEElement::MAX_NODES];
                for (int k=0; k<pme->Nodes(); ++k) cm[k] = m.Node(pme->m_node[k]).get(m_dofC + sid);
                double c2 = pme->eval(cm, data.m_rs[0], data.m_rs[1]);
                data.m_cg[m_ssl[isol]] = c1[isol] - c2;
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
                data.m_fls = m_bsmfls ? fls : ph/pn;
                data.m_mueff = m_mu*(1.0-(1.0-m_phi)*data.m_fls);
                data.m_mueff = MBRACKET(data.m_mueff);
                
                // total traction
                t = (nu + s1*data.m_mueff)*(-pn);
                
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
                    data.m_fls = m_bsmfls ? fls : ph/(-tn);
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
                        data.m_fls = m_bsmfls ? fls : ph/pn;
                    }
                    
                    // store the previous values as the current
                    data.m_pme = data.m_pmep;
                    data.m_rs = data.m_rsp;
                    
                    // recalculate gaps
                    data.m_dg = dg;
                    
                    // recalculate pressure gap
                    bool mporo = ms.m_bporo;
                    if (sporo && mporo)
                    {
                        double pm[FEElement::MAX_NODES];
                        for (int k=0; k<pme->Nodes(); ++k) pm[k] = m.Node(pme->m_node[k]).get(m_dofP);
                        double p2 = pme->eval(pm, data.m_rs[0], data.m_rs[1]);
                        data.m_pg = p - p2;
                    }
                    //TODO: what if stick/slip moves from poro to elastic or vice versa? Should I add
                    // an "else {data.m_pg = 0}" line?
                    
                    // recalculate concentration gaps
                    for (int isol=0; isol<nsol; ++isol) {
                        int sid = m_sid[isol];
                        double cm[FEElement::MAX_NODES];
                        for (int k=0; k<pme->Nodes(); ++k) cm[k] = m.Node(pme->m_node[k]).get(m_dofC + sid);
                        double c2 = pme->eval(cm, data.m_rs[0], data.m_rs[1]);
                        data.m_cg[m_ssl[isol]] = c1[isol] - c2;
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
                        data.m_fls = m_bsmfls ? fls : ph/pn;
                        data.m_mueff = m_mu*(1.0-(1.0-m_phi)*data.m_fls);
                        data.m_mueff = MBRACKET(data.m_mueff);
                        
                        // total traction
                        t = (nu + s1*data.m_mueff)*(-pn);
                        
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
                    data.m_fls = m_bsmfls ? fls : ph/pn;
                    data.m_mueff = m_mu*(1.0-(1.0-m_phi)*data.m_fls);
                    data.m_mueff = MBRACKET(data.m_mueff);
                    
                    // total traction
                    t = (nu + s1*data.m_mueff)*(-pn);
                    
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
void FESlidingInterfaceMP::LoadVector(FEGlobalVector& R, const FETimeInfo& tp)
{
    vector<int> sLM, mLM, LM, en;
    vector<double> fe;
    const int MN = FEElement::MAX_NODES;
    double detJ[MN], w[MN], *Hs, Hm[MN];
    double N[MN*10];
    int nsol = (int)m_sid.size();
    
    FEModel& fem = *GetFEModel();
    
    // need to multiply multiphasic stiffness entries by the timestep
    double dt = fem.GetTime().timeIncrement;

    double psf = GetPenaltyScaleFactor();
    
    m_ss.m_Ft = vec3d(0, 0, 0);
    m_ms.m_Ft = vec3d(0, 0, 0);
    
    // loop over the nr of passes
    int npass = (m_btwo_pass?2:1);
    for (int np=0; np<npass; ++np)
    {
        // get primary and secondary surface
        FESlidingSurfaceMP& ss = (np == 0? m_ss : m_ms);
        FESlidingSurfaceMP& ms = (np == 0? m_ms : m_ss);
        vector<int>& sl = (np == 0? m_ssl : m_msl);
        
        // loop over all primary surface elements
        for (int i=0; i<ss.Elements(); ++i)
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
                // get integration point data
                FESlidingSurfaceMP::Data& pt = static_cast<FESlidingSurfaceMP::Data&>(*se.GetMaterialPoint(j));
                
                // calculate contact pressure and account for stick
                double pn;
                vec3d t = ContactTraction(ss, i, j, ms, pn);
                
                // get the secondary surface element
                FESurfaceElement* pme = pt.m_pme;
                
                if (pme)
                {
                    // get the secondary surface element
                    FESurfaceElement& me = *pme;
                    
                    bool mporo = ms.m_bporo;
                    
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
                            // calculate nr of pressure dofs
                            int ndof = nseln + nmeln;
                            
                            // normal fluid flux
                            double epsp = m_epsp*pt.m_epsp*psf;
                            double wn = pt.m_Lmp + epsp*pt.m_pg;
                            
                            // fill the LM
                            LM.resize(ndof);
                            for (int k=0; k<nseln; ++k) LM[k      ] = sLM[3*nseln+k];
                            for (int k=0; k<nmeln; ++k) LM[k+nseln] = mLM[3*nmeln+k];
                            
                            // fill the force array
                            fe.resize(ndof);
                            zero(fe);
                            for (int k = 0; k<nseln; ++k) N[k      ] = Hs[k];
                            for (int k = 0; k<nmeln; ++k) N[k+nseln] = -Hm[k];
                            
                            for (int k = 0; k<ndof; ++k) fe[k] += dt*N[k]*wn*detJ[j]*w[j];
                            
                            
                            // assemble residual
                            R.Assemble(en, LM, fe);
                        }
                        
                        // do the solute stuff
                        for (int isol=0; isol<nsol; ++isol)
                        {
                            int sid = m_sid[isol];
                            
                            // calculate nr of concentration dofs
                            int ndof = nseln + nmeln;
                            
                            // calculate normal effective solute flux
                            int l = sl[isol];
                            double epsc = m_epsc*pt.m_epsc[l]*psf;
                            double jn = pt.m_Lmc[l] + epsc*pt.m_cg[l];
                            
                            // fill the LM
                            LM.resize(ndof);
                            for (int k = 0; k<nseln; ++k) LM[k      ] = sLM[(4+sid)*nseln+k];
                            for (int k = 0; k<nmeln; ++k) LM[k+nseln] = mLM[(4+sid)*nmeln+k];
                            
                            // fill the force array
                            fe.resize(ndof);
                            zero(fe);
                            for (int k = 0; k<nseln; ++k) N[k] = Hs[k];
                            for (int k = 0; k<nmeln; ++k) N[k+nseln] = -Hm[k];
                            
                            for (int k = 0; k<ndof; ++k) fe[k] += dt*N[k]*jn*detJ[j]*w[j];
                            
                            // assemble residual
                            R.Assemble(en, LM, fe);
                        }
                    }
                }
            }
        }
    }
}

//-----------------------------------------------------------------------------
void FESlidingInterfaceMP::StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp)
{
    int i, j, k, l;
    vector<int> sLM, mLM, LM, en;
    const int MN = FEElement::MAX_NODES;
    double detJ[MN], w[MN], *Hs, Hm[MN];
    FEElementMatrix ke;
    int nsol = (int)m_sid.size();
    vector<double> jn(nsol);
    
    FEModel& fem = *GetFEModel();
     
    double psf = GetPenaltyScaleFactor();
    
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
        FESlidingSurfaceMP& ss = (np == 0? m_ss : m_ms);
        FESlidingSurfaceMP& ms = (np == 0? m_ms : m_ss);
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
            
            double pn[MN] = {0};
            vector< vector<double> >cn(nsol,vector<double>(MN));
            if (sporo) {
                for (j=0; j<nseln; ++j)
                {
                    pn[j] = ss.GetMesh()->Node(se.m_node[j]).get(m_dofP);
                    for (int isol=0; isol<nsol; ++isol) {
                        cn[isol][j] = ss.GetMesh()->Node(se.m_node[j]).get(m_dofC + m_sid[isol]);
                    }
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
            }
            
            // loop over all integration points
            for (j=0; j<nint; ++j)
            {
                // get integration point data
                FESlidingSurfaceMP::Data& pt = static_cast<FESlidingSurfaceMP::Data&>(*se.GetMaterialPoint(j));
                
                // calculate contact traction and account for stick
                double pn;
                vec3d t = ContactTraction(ss, i, j, ms, pn);
                
                // get the secondary element
                FESurfaceElement* pme = pt.m_pme;
                
                // calculate normal effective solute flux
                for (int isol=0; isol<nsol; ++isol)
                {
                    int l = sl[isol];
                    double epsc = m_epsc*pt.m_epsc[l]*psf;
                    jn[isol] = pt.m_Lmc[l] + epsc*pt.m_cg[l];
                }
                
                // normal fluid flux
                double epsp = m_epsp*pt.m_epsp*psf;
                double wn = pt.m_Lmp + epsp*pt.m_pg;
                
                if (pme)
                {
                    FESurfaceElement& me = *pme;
                    
                    bool mporo = ms.m_bporo;
                    
                    // get the nr of secondary nodes
                    int nmeln = me.Nodes();
                    
                    // nodal data
                    double pm[MN] = {0};
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
                    
                    int ndpn;    // number of dofs per node
                    int ndof;    // number of dofs in stiffness matrix
                    
                    if (nsol) {
                        // calculate dofs for biphasic-solute contact
                        ndpn = 4+nsol;
                        ndof = ndpn*(nseln+nmeln);
                        
                        // build the LM vector
                        LM.resize(ndof);
                        
                        for (k=0; k<nseln; ++k)
                        {
                            LM[ndpn*k  ] = sLM[3*k  ];            // x-dof
                            LM[ndpn*k+1] = sLM[3*k+1];            // y-dof
                            LM[ndpn*k+2] = sLM[3*k+2];            // z-dof
                            LM[ndpn*k+3] = sLM[3*nseln+k];        // p-dof
                            for (int isol=0; isol<nsol; ++isol)
                                LM[ndpn*k+4+isol] = sLM[(4+m_sid[isol])*nseln+k];        // c-dof
                        }
                        for (k=0; k<nmeln; ++k)
                        {
                            LM[ndpn*(k+nseln)  ] = mLM[3*k  ];            // x-dof
                            LM[ndpn*(k+nseln)+1] = mLM[3*k+1];            // y-dof
                            LM[ndpn*(k+nseln)+2] = mLM[3*k+2];            // z-dof
                            LM[ndpn*(k+nseln)+3] = mLM[3*nmeln+k];        // p-dof
                            for (int isol=0; isol<nsol; ++isol)
                                LM[ndpn*(k+nseln)+4+isol] = mLM[(4+m_sid[isol])*nmeln+k];        // c-dof
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
                            LM[ndpn*k  ] = sLM[3*k  ];            // x-dof
                            LM[ndpn*k+1] = sLM[3*k+1];            // y-dof
                            LM[ndpn*k+2] = sLM[3*k+2];            // z-dof
                            LM[ndpn*k+3] = sLM[3*nseln+k];        // p-dof
                        }
                        for (k=0; k<nmeln; ++k)
                        {
                            LM[ndpn*(k+nseln)  ] = mLM[3*k  ];            // x-dof
                            LM[ndpn*(k+nseln)+1] = mLM[3*k+1];            // y-dof
                            LM[ndpn*(k+nseln)+2] = mLM[3*k+2];            // z-dof
                            LM[ndpn*(k+nseln)+3] = mLM[3*nmeln+k];        // p-dof
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
                    
                    // primary shape functions
                    Hs = se.H(j);
                    
                    // secondary shape functions
                    double r = pt.m_rs[0];
                    double s = pt.m_rs[1];
                    me.shape_fnc(Hm, r, s);
                    
                    // get primary normal vector
                    vec3d nu = pt.m_nu;
                    
                    // gap function
                    double g = pt.m_gap;
                    
                    // penalty
                    double eps = m_epsn*pt.m_epsn*psf;
                    
                    // only evaluate stiffness matrix if contact traction is non-zero
                    if (pn > 0) {
                        // if stick
                        if (pt.m_bstick) {
                            
                            // create the stiffness matrix
                            ke.resize(ndof, ndof); ke.zero();
                            
                            // evaluate basis vectors on primary surface
                            vec3d gscov[2];
                            ss.CoBaseVectors(se, j, gscov);
                            
                            // identity tensor
                            mat3d I = mat3dd(1);
                            
                            // evaluate Mc and Ac and combine them into As
                            double* Gr = se.Gr(j);
                            double* Gs = se.Gs(j);
                            mat3d Ac[MN], As[MN];
                            mat3d gscovh[2];
                            gscovh[0].skew(gscov[0]); gscovh[1].skew(gscov[1]);
                            for (int k=0; k<nseln; ++k) {
                                Ac[k] = (gscovh[1]*Gr[k] - gscovh[0]*Gs[k])/detJ[j];
                                As[k] = t & (Ac[k]*nu);
                            }
                            
                            // --- S O L I D - S O L I D   C O N T A C T ---
                            
                            double tmp = detJ[j]*w[j];
                            for (int a=0; a<nseln; ++a) {
                                k = a*ndpn;
                                for (int c=0; c<nseln; ++c) {
                                    l = c*ndpn;
                                    mat3d Kac = (I*Hs[a]*Hs[c]*eps - As[c]*Hs[a])*tmp;
                                    ke[k  ][l  ] += Kac(0,0); ke[k  ][l+1] += Kac(0,1); ke[k  ][l+2] += Kac(0,2);
                                    ke[k+1][l  ] += Kac(1,0); ke[k+1][l+1] += Kac(1,1); ke[k+1][l+2] += Kac(1,2);
                                    ke[k+2][l  ] += Kac(2,0); ke[k+2][l+1] += Kac(2,1); ke[k+2][l+2] += Kac(2,2);
                                }
                                for (int d=0; d<nmeln; ++d) {
                                    l = (nseln+d)*ndpn;
                                    mat3d Kad = I*(-Hs[a]*Hm[d]*eps)*tmp;
                                    ke[k  ][l  ] += Kad(0,0); ke[k  ][l+1] += Kad(0,1); ke[k  ][l+2] += Kad(0,2);
                                    ke[k+1][l  ] += Kad(1,0); ke[k+1][l+1] += Kad(1,1); ke[k+1][l+2] += Kad(1,2);
                                    ke[k+2][l  ] += Kad(2,0); ke[k+2][l+1] += Kad(2,1); ke[k+2][l+2] += Kad(2,2);
                                }
                            }
                            for (int b=0; b<nmeln; ++b) {
                                k = (nseln+b)*ndpn;
                                for (int c=0; c<nseln; ++c) {
                                    l = c*ndpn;
                                    mat3d Kbc = (I*(-Hm[b]*Hs[c]*eps) + As[c]*Hm[b])*tmp;
                                    ke[k  ][l  ] += Kbc(0,0); ke[k  ][l+1] += Kbc(0,1); ke[k  ][l+2] += Kbc(0,2);
                                    ke[k+1][l  ] += Kbc(1,0); ke[k+1][l+1] += Kbc(1,1); ke[k+1][l+2] += Kbc(1,2);
                                    ke[k+2][l  ] += Kbc(2,0); ke[k+2][l+1] += Kbc(2,1); ke[k+2][l+2] += Kbc(2,2);
                                }
                                for (int d=0; d<nmeln; ++d) {
                                    l = (nseln+d)*ndpn;
                                    mat3d Kbd = I*Hm[b]*Hm[d]*eps*tmp;
                                    ke[k  ][l  ] += Kbd(0,0); ke[k  ][l+1] += Kbd(0,1); ke[k  ][l+2] += Kbd(0,2);
                                    ke[k+1][l  ] += Kbd(1,0); ke[k+1][l+1] += Kbd(1,1); ke[k+1][l+2] += Kbd(1,2);
                                    ke[k+2][l  ] += Kbd(2,0); ke[k+2][l+1] += Kbd(2,1); ke[k+2][l+2] += Kbd(2,2);
                                }
                            }
                            
                            // --- M U L T I P H A S I C   S T I F F N E S S ---
                            if (sporo && mporo)
                            {
                                // need to multiply biphasic stiffness entries by the timestep
                                double dt = fem.GetTime().timeIncrement;
                                
                                // --- S O L I D - P R E S S U R E / S O L U T E   C O N T A C T ---
                                
                                for (int a=0; a<nseln; ++a) {
                                    k = a*ndpn;
                                    for (int c=0; c<nseln; ++c) {
                                        l = c*ndpn;
                                        vec3d gac = (Ac[c]*nu)*(-Hs[a]*wn)*tmp*dt;
                                        ke[k+3][l  ] += gac.x; ke[k+3][l+1] += gac.y; ke[k+3][l+2] += gac.z;
                                        for (int isol=0; isol<nsol; ++isol) {
                                            vec3d hac = (Ac[c]*nu)*(-Hs[a]*jn[isol])*tmp*dt;
                                            ke[k+4+isol][l  ] += hac.x; ke[k+4+isol][l+1] += hac.y; ke[k+4+isol][l+2] += hac.z;
                                        }
                                    }
                                    for (int d=0; d<nmeln; ++d) {
                                        l = (nseln+d)*ndpn;
                                        vec3d gad = vec3d(0,0,0);
                                        ke[k+3][l  ] += gad.x; ke[k+3][l+1] += gad.y; ke[k+3][l+2] += gad.z;
                                        for (int isol=0; isol<nsol; ++isol) {
                                            vec3d had = vec3d(0,0,0);
                                            ke[k+4+isol][l  ] += had.x; ke[k+4+isol][l+1] += had.y; ke[k+4+isol][l+2] += had.z;
                                        }
                                    }
                                }
                                for (int b=0; b<nmeln; ++b) {
                                    k = (nseln+b)*ndpn;
                                    for (int c=0; c<nseln; ++c) {
                                        l = c*ndpn;
                                        vec3d gbc = (Ac[c]*nu)*Hm[b]*wn*tmp*dt;
                                        ke[k+3][l  ] += gbc.x; ke[k+3][l+1] += gbc.y; ke[k+3][l+2] += gbc.z;
                                        for (int isol=0; isol<nsol; ++isol) {
                                            vec3d hbc = (Ac[c]*nu)*Hm[b]*jn[isol]*tmp*dt;
                                            ke[k+4+isol][l  ] += hbc.x; ke[k+4+isol][l+1] += hbc.y; ke[k+4+isol][l+2] += hbc.z;
                                        }
                                    }
                                    for (int d=0; d<nmeln; ++d) {
                                        l = (nseln+d)*ndpn;
                                        vec3d gbd = vec3d(0,0,0);
                                        ke[k+3][l  ] += gbd.x; ke[k+3][l+1] += gbd.y; ke[k+3][l+2] += gbd.z;
                                        for (int isol=0; isol<nsol; ++isol) {
                                            vec3d hbd = vec3d(0,0,0);
                                            ke[k+4+isol][l  ] += hbd.x; ke[k+4+isol][l+1] += hbd.y; ke[k+4+isol][l+2] += hbd.z;
                                        }
                                    }
                                }
                                
                                // --- P R E S S U R E - P R E S S U R E  /  C O N C E N T R A T I O N - C O N C E N T R A T I O N  C O N T A C T ---
                                
                                for (int a=0; a<nseln; ++a) {
                                    k = a*ndpn;
                                    for (int c=0; c<nseln; ++c) {
                                        l = c*ndpn;
                                        double gac = (-Hs[a]*Hs[c]*epsp)*tmp*dt;
                                        ke[k+3][l+3] += gac;
                                        for (int isol=0; isol<nsol; ++isol) {
                                            for (int jsol=0; jsol<nsol; ++jsol) {
                                                int z = (isol == jsol? 1.0 : 0.0);
                                                double epsc = m_epsc*pt.m_epsc[sl[isol]]*psf;
                                                double hac = (-Hs[a]*Hs[c]*epsc*z)*tmp*dt;
                                                ke[k+4+isol][l+4+jsol] += hac;
                                            }
                                        }
                                    }
                                    for (int d=0; d<nmeln; ++d) {
                                        l = (nseln+d)*ndpn;
                                        double gad = (Hs[a]*Hm[d]*epsp)*tmp*dt;
                                        ke[k+3][l+3] += gad;
                                        for (int isol=0; isol<nsol; ++isol) {
                                            for (int jsol=0; jsol<nsol; ++jsol) {
                                                int z = (isol == jsol? 1.0 : 0.0);
                                                double epsc = m_epsc*pt.m_epsc[sl[isol]]*psf;
                                                double had = (Hs[a]*Hm[d]*epsc*z)*tmp*dt;
                                                ke[k+4+isol][l+4+jsol] += had;
                                            }
                                        }
                                    }
                                }
                                for (int b=0; b<nmeln; ++b) {
                                    k = (nseln+b)*ndpn;
                                    for (int c=0; c<nseln; ++c) {
                                        l = c*ndpn;
                                        double gbc = (Hm[b]*Hs[c]*epsp)*tmp*dt;
                                        ke[k+3][l+3] += gbc;
                                        for (int isol=0; isol<nsol; ++isol) {
                                            for (int jsol=0; jsol<nsol; ++jsol) {
                                                int z = (isol == jsol? 1.0 : 0.0);
                                                double epsc = m_epsc*pt.m_epsc[sl[isol]]*psf;
                                                double hbc = (Hm[b]*Hs[c]*epsc*z)*tmp*dt;
                                                ke[k+4+isol][l+4+jsol] += hbc;
                                            }
                                        }
                                    }
                                    for (int d=0; d<nmeln; ++d) {
                                        l = (nseln+d)*ndpn;
                                        double gbd = (-Hm[b]*Hm[d]*epsp)*tmp*dt;
                                        ke[k+3][l+3] += gbd;
                                        for (int isol=0; isol<nsol; ++isol) {
                                            for (int jsol=0; jsol<nsol; ++jsol) {
                                                int z = (isol == jsol? 1.0 : 0.0);
                                                double epsc = m_epsc*pt.m_epsc[sl[isol]]*psf;
                                                double hbd = (-Hm[b]*Hm[d]*epsc*z)*tmp*dt;
                                                ke[k+4+isol][l+4+jsol] += hbd;
                                            }
                                        }
                                    }
                                }
                            }
                            
                            // assemble the global stiffness
                            ke.SetNodes(en);
                            ke.SetIndices(LM);
                            LS.Assemble(ke);
                        }
                        // if slip
                        else {
                            
                            // create the stiffness matrix
                            ke.resize(ndof, ndof); ke.zero();
                            
                            double tn = -pn;
                            
                            // obtain the slip direction s1 and inverse of spatial increment dh
                            double dh = 0, hd = 0;
                            vec3d dr(0,0,0);
                            vec3d s1 = SlipTangent(ss, i, j, ms, dh, dr);
                            
                            if (dh != 0) hd = 1.0 / dh;
                            
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
                            mat3d Pn = mat3dd(1) - (nu & nu);
                            mat3d Nb1 = mat3dd(1) - (gscov[0] & gscnt[0]) - (gscov[1] & gscnt[1]);
                            mat3d Nt1 = nu & (Nb1*nu);
                            mat3d S1 = s1 & nu;
                            mat3d Ps = (mat3dd(1) - (s1 & s1))*hd;
                            mat3d St1 = s1 & (Nb1*nu);
                            
                            // evaluate frictional contact vectors and tensors
                            vec3d m = ((dgscov[0] ^ gscov[1]) + (gscov[0] ^ dgscov[1]));
                            vec3d c1 = Pn*m*(1/detJ[j]);
                            mat3d Q1 = (mat3dd(1)*(nu * m) + (nu & m))*(1/detJ[j]);
                            mat3d B = ((Ps*c1) & (Nb1*nu)) - Ps*Pn;
                            mat3d R = (mat3dd(1)*(nu * dr) + (nu & dr))/(-g);
                            mat3d L1 = Ps*(Pn*Q1 + R - mat3dd(1))*Pn*(-g);
                            
                            // evaluate Ac, Mc, and combine into As
                            double* Gr = se.Gr(j);
                            double* Gs = se.Gs(j);
                            mat3d gscovh[2];
                            mat3d dgscovh[2];
                            gscovh[0].skew(gscov[0]); gscovh[1].skew(gscov[1]);
                            dgscovh[0].skew(dgscov[0]); dgscovh[1].skew(dgscov[1]);
                            mat3d Ac[MN];
                            mat3d As[MN];
                            mat3d Pc[MN];
                            for (int c=0; c<nseln; ++c){
                                vec3d mc = gscnt[0]*Gr[c] + gscnt[1]*Gs[c];
                                mat3d Mc = nu & mc;
                                Ac[c] = (gscovh[1]*Gr[c] - gscovh[0]*Gs[c])/detJ[j];
                                mat3d Acb = (dgscovh[1]*Gr[c] - dgscovh[0]*Gs[c])/detJ[j];
                                vec3d hcp = (N1*mc + Ac[c]*nu);
                                vec3d hcmb = (N1*mc*m_mu - Ac[c]*nu*pt.m_mueff);
                                As[c] = Ac[c]+Mc*N1;
                                mat3d Jc = (L1*Ac[c]) - Ps*Pn*Acb*(-g);
                                Pc[c] = (s1 & hcmb) + ((Ps*c1) & hcp)*pt.m_mueff*(-g) - Jc*pt.m_mueff;
                            }
                            
                            // evaluate mb and Mb
                            double Hmr[MN], Hms[MN];
                            me.shape_deriv(Hmr, Hms, r, s);
                            vec3d mb[MN];
                            mat3d Pb[MN];
                            for (k=0; k<nmeln; ++k) {
                                mb[k] = gmcnt[0]*Hmr[k] + gmcnt[1]*Hms[k];
                                Pb[k] = ((-nu) & mb[k]) - (s1 & mb[k])*pt.m_mueff;
                            }
                            
                            // evaluate Gbc
                            matrix Gbc(nmeln,nseln);
                            for (int b=0; b<nmeln; ++b) {
                                for (int c=0; c<nseln; ++c) {
                                    Gbc(b,c)
                                    = (a[0][0]*Hmr[b]*Gr[c]
                                       + a[0][1]*Hmr[b]*Gs[c]
                                       + a[1][0]*Hms[b]*Gr[c]
                                       + a[1][1]*Hms[b]*Gs[c])*(-g);
                                }
                            }
                            
                            // define Tt, T
                            mat3d Tt = Nt1 + St1*m_mu;
                            mat3d T = N1 + S1*pt.m_mueff;
                            
                            // --- S O L I D - S O L I D   C O N T A C T ---
                            
                            // All Kmn terms have the opposite sign from Brandon's notes
                            // Does FEBio put the negative sign in the stiffness matrix? K*u = f => f - K*u = 0?
                            
                            double tmp = detJ[j]*w[j];
                            for (int a=0; a<nseln; ++a) {
                                k = a*ndpn;
                                for (int c=0; c<nseln; ++c) {
                                    l = c*ndpn;
                                    mat3d Kac = ((Tt*eps+B*tn*pt.m_mueff)*Hs[a]*Hs[c]+(As[c]+Pc[c])*Hs[a]*tn)*tmp;
                                    ke[k  ][l  ] += Kac(0,0); ke[k  ][l+1] += Kac(0,1); ke[k  ][l+2] += Kac(0,2);
                                    ke[k+1][l  ] += Kac(1,0); ke[k+1][l+1] += Kac(1,1); ke[k+1][l+2] += Kac(1,2);
                                    ke[k+2][l  ] += Kac(2,0); ke[k+2][l+1] += Kac(2,1); ke[k+2][l+2] += Kac(2,2);
                                }
                                for (int d=0; d<nmeln; ++d) {
                                    l = (nseln+d)*ndpn;
                                    mat3d Kad = (Tt*eps+B*tn*pt.m_mueff)*(-Hs[a]*Hm[d]*tmp);
                                    ke[k  ][l  ] += Kad(0,0); ke[k  ][l+1] += Kad(0,1); ke[k  ][l+2] += Kad(0,2);
                                    ke[k+1][l  ] += Kad(1,0); ke[k+1][l+1] += Kad(1,1); ke[k+1][l+2] += Kad(1,2);
                                    ke[k+2][l  ] += Kad(2,0); ke[k+2][l+1] += Kad(2,1); ke[k+2][l+2] += Kad(2,2);
                                }
                            }
                            for (int b=0; b<nmeln; ++b) {
                                k = (nseln+b)*ndpn;
                                for (int c=0; c<nseln; ++c) {
                                    l = c*ndpn;
                                    mat3d Kbc = ((Tt*eps+B*tn*pt.m_mueff)*(-Hm[b]*Hs[c])-(As[c]+Pc[c])*(Hm[b]*tn)-Pb[b]*(Hs[c]*tn)-T*(Gbc[b][c]*tn))*tmp;
                                    ke[k  ][l  ] += Kbc(0,0); ke[k  ][l+1] += Kbc(0,1); ke[k  ][l+2] += Kbc(0,2);
                                    ke[k+1][l  ] += Kbc(1,0); ke[k+1][l+1] += Kbc(1,1); ke[k+1][l+2] += Kbc(1,2);
                                    ke[k+2][l  ] += Kbc(2,0); ke[k+2][l+1] += Kbc(2,1); ke[k+2][l+2] += Kbc(2,2);
                                }
                                for (int d=0; d<nmeln; ++d) {
                                    l = (nseln+d)*ndpn;
                                    mat3d Kbd = ((Tt*eps+B*tn*pt.m_mueff)*Hm[b]*Hm[d]+Pb[b]*Hm[d]*tn)*tmp;
                                    ke[k  ][l  ] += Kbd(0,0); ke[k  ][l+1] += Kbd(0,1); ke[k  ][l+2] += Kbd(0,2);
                                    ke[k+1][l  ] += Kbd(1,0); ke[k+1][l+1] += Kbd(1,1); ke[k+1][l+2] += Kbd(1,2);
                                    ke[k+2][l  ] += Kbd(2,0); ke[k+2][l+1] += Kbd(2,1); ke[k+2][l+2] += Kbd(2,2);
                                }
                            }
                            
                            // --- M U L T I P H A S I C   S T I F F N E S S ---
                            if (sporo && mporo)
                            {
                                double dt = fem.GetTime().timeIncrement;
                                
                                double epsp = m_epsp*pt.m_epsp*psf;
                                
                                // p vector (gradients of effective solute concentrations on master surface)
                                double dpmr = me.eval_deriv1(pm, r, s);
                                double dpms = me.eval_deriv2(pm, r, s);
                                vec3d p = gmcnt[0]*dpmr + gmcnt[1]*dpms;
                                
                                // evaluate Pc
                                double Pc[MN];
                                for (int k=0; k<nseln; ++k) {
                                    Pc[k] = (a[0][0]*dpmr*Gr[k]
                                             + a[0][1]*dpmr*Gs[k]
                                             + a[1][0]*dpms*Gr[k]
                                             + a[1][1]*dpms*Gs[k])*(-g);
                                }
                                
                                // q vectors (gradients of effective solute concentrations on master surface)
                                // Cc scalars
                                vector<vec3d> q(nsol);
                                vector< vector<double> > Cc(nsol, vector<double>(MN));
                                for (int isol=0; isol<nsol; ++isol) {
                                    double dcmr = me.eval_deriv1(&cm[isol][0], r, s);
                                    double dcms = me.eval_deriv2(&cm[isol][0], r, s);
                                    q[isol] = gmcnt[0]*dcmr + gmcnt[1]*dcms;
                                    for (int k=0; k<nseln; ++k) {
                                        Cc[isol][k] = (a[0][0]*dcmr*Gr[k]
                                                       + a[0][1]*dcmr*Gs[k]
                                                       + a[1][0]*dcms*Gr[k]
                                                       + a[1][1]*dcms*Gs[k])*(-g);
                                    }
                                }
                                
                                // --- S O L I D - P R E S S U R E / S O L U T E   C O N T A C T ---
                                
                                for (int a=0; a<nseln; ++a) {
                                    k = a*ndpn;
                                    for (int c=0; c<nseln; ++c) {
                                        l = c*ndpn;
                                        vec3d gac = (p*Hs[a]*Hs[c]*epsp-((Ac[c]*nu)*wn+nu*epsp*Pc[c])*Hs[a])*tmp*dt;
                                        ke[k+3][l  ] += gac.x; ke[k+3][l+1] += gac.y; ke[k+3][l+2] += gac.z;
                                        vec3d kac = (s1*m_mu*(1.0-m_phi))*(-Hs[a]*Hs[c])*tmp*dt;
                                        ke[k  ][l+3] += kac.x; ke[k+1][l+3] += kac.y; ke[k+2][l+3] += kac.z;
                                        for (int isol=0; isol<nsol; ++isol) {
                                            double epsc = m_epsc*pt.m_epsc[sl[isol]]*psf;
                                            vec3d hac = (q[isol]*Hs[a]*Hs[c]*epsc-((Ac[c]*nu)*jn[isol]+nu*Cc[isol][c]*epsc)*Hs[a])*tmp*dt;
                                            ke[k+4+isol][l  ] += hac.x; ke[k+4+isol][l+1] += hac.y; ke[k+4+isol][l+2] += hac.z;
                                        }
                                    }
                                    for (int d=0; d<nmeln; ++d) {
                                        l = (nseln+d)*ndpn;
                                        vec3d gad = p*(-Hs[a]*Hm[d]*epsp)*tmp*dt;
                                        ke[k+3][l  ] += gad.x; ke[k+3][l+1] += gad.y; ke[k+3][l+2] += gad.z;
                                        for (int isol=0; isol<nsol; ++isol) {
                                            double epsc = m_epsc*pt.m_epsc[sl[isol]]*psf;
                                            vec3d had = q[isol]*(-Hs[a]*Hm[d]*epsc)*tmp*dt;
                                            ke[k+4+isol][l  ] += had.x; ke[k+4+isol][l+1] += had.y; ke[k+4+isol][l+2] += had.z;
                                        }
                                    }
                                }
                                for (int b=0; b<nmeln; ++b) {
                                    k = (nseln+b)*ndpn;
                                    for (int c=0; c<nseln; ++c) {
                                        l = c*ndpn;
                                        vec3d gbc = (p*(-Hm[b]*Hs[c]*epsp)+((Ac[c]*nu)*wn + nu*Pc[c]*epsp)*Hm[b] + (mb[b]*wn*Hs[c]) - nu*Gbc[b][c]*wn)*tmp*dt;
                                        ke[k+3][l  ] += gbc.x; ke[k+3][l+1] += gbc.y; ke[k+3][l+2] += gbc.z;
                                        vec3d kbc = (s1*m_mu*(1.0-m_phi))*(Hm[b]*Hs[c])*tmp*dt;
                                        ke[k  ][l+3] += kbc.x; ke[k+1][l+3] += kbc.y; ke[k+2][l+3] += kbc.z;
                                        for (int isol=0; isol<nsol; ++isol) {
                                            double epsc = m_epsc*pt.m_epsc[sl[isol]]*psf;
                                            vec3d hbc = (q[isol]*(-Hm[b]*Hs[c]*epsc)+((Ac[c]*nu)*jn[isol]+nu*epsc*Cc[isol][c])*Hm[b]+(mb[b]*jn[isol]*Hs[c])-(nu*Gbc[b][c]*jn[isol]))*tmp*dt;
                                            ke[k+4+isol][l  ] += hbc.x; ke[k+4+isol][l+1] += hbc.y; ke[k+4+isol][l+2] += hbc.z;
                                        }
                                    }
                                    for (int d=0; d<nmeln; ++d) {
                                        l = (nseln+d)*ndpn;
                                        vec3d gbd = (p*(Hm[b]*Hm[d]*epsp)-(mb[b]*wn*Hm[d]))*tmp*dt;
                                        ke[k+3][l  ] += gbd.x; ke[k+3][l+1] += gbd.y; ke[k+3][l+2] += gbd.z;
                                        for (int isol=0; isol<nsol; ++isol) {
                                            double epsc = m_epsc*pt.m_epsc[sl[isol]]*psf;
                                            vec3d hbd = (q[isol]*(Hm[b]*Hm[d]*epsc)-(mb[b]*jn[isol]*Hm[d]))*tmp*dt;
                                            ke[k+4+isol][l  ] += hbd.x; ke[k+4+isol][l+1] += hbd.y; ke[k+4+isol][l+2] += hbd.z;
                                        }
                                    }
                                }
                                
                                // --- P R E S S U R E - P R E S S U R E  /  C O N C E N T R A T I O N - C O N C E N T R A T I O N  C O N T A C T ---
                                
                                for (int a=0; a<nseln; ++a) {
                                    k = a*ndpn;
                                    for (int c=0; c<nseln; ++c) {
                                        l = c*ndpn;
                                        double gac = (-Hs[a]*Hs[c]*epsp)*tmp*dt;
                                        ke[k+3][l+3] += gac;
                                        for (int isol=0; isol<nsol; ++isol) {
                                            for (int jsol=0; jsol<nsol; ++jsol) {
                                                int z = (isol == jsol? 1.0 : 0.0);
                                                double epsc = m_epsc*pt.m_epsc[sl[isol]]*psf;
                                                double hac = (-Hs[a]*Hs[c]*epsc*z)*(dt*detJ[j]*w[j]);
                                                ke[k+4+isol][l+4+jsol] += hac;
                                            }
                                        }
                                    }
                                    for (int d=0; d<nmeln; ++d) {
                                        l = (nseln+d)*ndpn;
                                        double gad = (Hs[a]*Hm[d]*epsp)*(dt*detJ[j]*w[j]);
                                        ke[k+3][l+3] += gad;
                                        for (int isol=0; isol<nsol; ++isol) {
                                            for (int jsol=0; jsol<nsol; ++jsol) {
                                                int z = (isol == jsol? 1.0 : 0.0);
                                                double epsc = m_epsc*pt.m_epsc[sl[isol]]*psf;
                                                double had = (Hs[a]*Hm[d]*epsc*z)*(dt*detJ[j]*w[j]);
                                                ke[k+4+isol][l+4+jsol] += had;
                                            }
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
                                            for (int jsol=0; jsol<nsol; ++jsol) {
                                                int z = (isol == jsol? 1.0 : 0.0);
                                                double epsc = m_epsc*pt.m_epsc[sl[isol]]*psf;
                                                double hbc = (Hm[b]*Hs[c]*epsc*z)*(dt*detJ[j]*w[j]);
                                                ke[k+4+isol][l+4+jsol] += hbc;
                                            }
                                        }
                                    }
                                    for (int d=0; d<nmeln; ++d) {
                                        l = (nseln+d)*ndpn;
                                        double gbd = (-Hm[b]*Hm[d]*epsp)*(dt*detJ[j]*w[j]);
                                        ke[k+3][l+3] += gbd;
                                        for (int isol=0; isol<nsol; ++isol) {
                                            for (int jsol=0; jsol<nsol; ++jsol) {
                                                int z = (isol == jsol? 1.0 : 0.0);
                                                double epsc = m_epsc*pt.m_epsc[sl[isol]]*psf;
                                                double hbd = (-Hm[b]*Hm[d]*epsc*z)*(dt*detJ[j]*w[j]);
                                                ke[k+4+isol][l+4+jsol] += hbd;
                                            }
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
    }
}

//-----------------------------------------------------------------------------
void FESlidingInterfaceMP::UpdateContactPressures()
{
    int np, n, i, j;
    const int MN = FEElement::MAX_NODES;
    const int MI = FEElement::MAX_INTPOINTS;
    double psf = GetPenaltyScaleFactor();
    int npass = (m_btwo_pass?2:1);
    for (np=0; np<npass; ++np)
    {
        FESlidingSurfaceMP& ss = (np == 0? m_ss : m_ms);
        FESlidingSurfaceMP& ms = (np == 0? m_ms : m_ss);
        
        // loop over all elements of the primary surface
        for (n=0; n<ss.Elements(); ++n)
        {
            FESurfaceElement& el = ss.Element(n);
            int nint = el.GaussPoints();
            
            // get the normal tractions at the integration points
            for (i=0; i<nint; ++i)
            {
                // get integration point data
                FESlidingSurfaceMP::Data& sd = static_cast<FESlidingSurfaceMP::Data&>(*el.GetMaterialPoint(i));
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
                    int mint = pme->GaussPoints();
                    vec3d ti[MI];
                    double pi[MI];
                    for (j=0; j<mint; ++j) {
                        FESlidingSurfaceMP::Data& md = static_cast<FESlidingSurfaceMP::Data&>(*pme->GetMaterialPoint(j));
                        
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
                            ti[j] = -(md.m_nu + md.m_s1*md.m_mueff)*pi[j];
                        }
                    }
                    // project the data to the nodes
                    vec3d tn[MN];
                    double pn[MN];
                    pme->FEElement::project_to_nodes(pi, pn);
                    pme->project_to_nodes(ti, tn);
                    // now evaluate the traction at the intersection point
                    double Ln = pme->eval(pn, sd.m_rs[0], sd.m_rs[1]);
                    vec3d trac = pme->eval(tn, sd.m_rs[0], sd.m_rs[1]);
                    sd.m_Ln += MBRACKET(Ln);
                    // tractions on secondary-primary are opposite, so subtract
                    sd.m_tr -= trac;
                }
            }
        }
        ss.EvaluateNodalContactPressures();
        ss.EvaluateNodalContactTractions();
    }
}

//-----------------------------------------------------------------------------
bool FESlidingInterfaceMP::Augment(int naug, const FETimeInfo& tp)
{
    // make sure we need to augment
    if (m_laugon != 1) return true;

    double Ln, Lp;
    int nsol = (int)m_sid.size();
    vector<double>Lc(nsol);
    bool bconv = true;
    
    double psf = GetPenaltyScaleFactor();
    
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
        for (int j=0; j<el.GaussPoints(); ++j)
        {
            FESlidingSurfaceMP::Data& ds = static_cast<FESlidingSurfaceMP::Data&>(*el.GetMaterialPoint(j));
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
            FESlidingSurfaceMP::Data& dm = static_cast<FESlidingSurfaceMP::Data&>(*el.GetMaterialPoint(j));
            if (dm.m_bstick)
                normL0 += dm.m_Lmt*dm.m_Lmt;
            else
                normL0 += dm.m_Lmd*dm.m_Lmd;
        }
    }
    
    // b. gap component
    // (is calculated during update)
    double maxgap = 0, maxpg = 0;
    vector<double> maxcg(nsol,0);
    
    // update Lagrange multipliers
    double normL1 = 0, eps, epsp, epsc;
    for (int i=0; i<m_ss.Elements(); ++i) {
        FESurfaceElement& el = m_ss.Element(i);
        vec3d tn[FEElement::MAX_INTPOINTS];
        if (m_bsmaug) m_ss.GetGPSurfaceTraction(i, tn);
        for (int j=0; j<el.GaussPoints(); ++j)
        {
            FESlidingSurfaceMP::Data& ds = static_cast<FESlidingSurfaceMP::Data&>(*el.GetMaterialPoint(j));
            
            // update Lagrange multipliers on primary surface
            eps = m_epsn*ds.m_epsn*psf;
            if (ds.m_bstick) {
                // if stick, augment total traction
                if (m_bsmaug) {
                    // replace this multiplier with a smoother version
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
                    // replace this multiplier with a smoother version
                    Ln = -(tn[j]*ds.m_nu);
                    ds.m_Lmd = MBRACKET(Ln);
                    if (m_btwo_pass) ds.m_Lmd /= 2;
                }
                else {
                    Ln = ds.m_Lmd + eps*ds.m_gap;
                    ds.m_Lmd = MBRACKET(Ln);
                }
                // then derive total traction
                double mueff = m_mu*(1.0-(1.0-m_phi)*(ds.m_p1-m_ambp)/ds.m_Lmd);
                mueff = MBRACKET(mueff);
                ds.m_Lmt = -(ds.m_nu*ds.m_Lmd + ds.m_s1*ds.m_Lmd*mueff);
                normL1 += ds.m_Lmd*ds.m_Lmd;
                
                if (Ln > 0) maxgap = max(maxgap, fabs(ds.m_gap));
            }
            
            if (m_ss.m_bporo) {
                Lp = 0;
                Lc.assign(nsol, 0);
                if (Ln > 0) {
                    epsp = m_epsp*ds.m_epsp*psf;
                    Lp = ds.m_Lmp + epsp*ds.m_pg;
                    maxpg = max(maxpg,fabs(ds.m_pg));
                    normDP += ds.m_pg*ds.m_pg;
                    for (int isol=0; isol<nsol; ++isol) {
                        int l = m_ssl[isol];
                        epsc = m_epsc*ds.m_epsc[l]*psf;
                        Lc[isol] = ds.m_Lmc[l] + epsc*ds.m_cg[l];
                        maxcg[isol] = max(maxcg[isol],fabs(ds.m_cg[l]));
                        normDC[isol] += ds.m_cg[l]*ds.m_cg[l];
                    }
                }
                ds.m_Lmp = Lp;
                for (int isol=0; isol<nsol; ++isol) ds.m_Lmc[m_ssl[isol]] = Lc[isol];
            }
        }
    }
    
    for (int i=0; i<m_ms.Elements(); ++i) {
        FESurfaceElement& el = m_ms.Element(i);
        vec3d tn[FEElement::MAX_INTPOINTS];
        if (m_bsmaug) m_ms.GetGPSurfaceTraction(i, tn);
        for (int j=0; j<el.GaussPoints(); ++j) {
            FESlidingSurfaceMP::Data& dm = static_cast<FESlidingSurfaceMP::Data&>(*el.GetMaterialPoint(j));
            
            // update Lagrange multipliers on master surface
            double eps = m_epsn*dm.m_epsn*psf;
            if (dm.m_bstick) {
                // if stick, augment total traction
                if (m_bsmaug) {
                    // replace this multiplier with a smoother version
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
                    // replace this multiplier with a smoother version
                    Ln = -(tn[j]*dm.m_nu);
                    dm.m_Lmd = MBRACKET(Ln);
                    if (m_btwo_pass) dm.m_Lmd /= 2;
                }
                else {
                    Ln = dm.m_Lmd + eps*dm.m_gap;
                    dm.m_Lmd = MBRACKET(Ln);
                }
                // then derive total traction
                double mueff = m_mu*(1.0-(1.0-m_phi)*(dm.m_p1-m_ambp)/dm.m_Lmd);
                mueff = MBRACKET(mueff);
                dm.m_Lmt = -(dm.m_nu*dm.m_Lmd + dm.m_s1*dm.m_Lmd*mueff);
                normL1 += dm.m_Lmd*dm.m_Lmd;
                
                if (Ln > 0) maxgap = max(maxgap, fabs(dm.m_gap));
            }
            
            if (m_ms.m_bporo) {
                Lp = 0;
                Lc.assign(nsol, 0);
                if (Ln > 0) {
                    epsp = m_epsp*dm.m_epsp*psf;
                    Lp = dm.m_Lmp + epsp*dm.m_pg;
                    maxpg = max(maxpg,fabs(dm.m_pg));
                    normDP += dm.m_pg*dm.m_pg;
                    for (int isol=0; isol<nsol; ++isol) {
                        int l = m_ssl[isol];
                        epsc = m_epsc*dm.m_epsc[l]*psf;
                        Lc[isol] = dm.m_Lmc[l] + epsc*dm.m_cg[l];
                        maxcg[isol] = max(maxcg[isol],fabs(dm.m_cg[l]));
                        normDC[isol] += dm.m_cg[l]*dm.m_cg[l];
                    }
                }
                dm.m_Lmp = Lp;
                for (int isol=0; isol<nsol; ++isol) dm.m_Lmc[m_ssl[isol]] = Lc[isol];
            }
        }
    }
    
    // normP should be a measure of the fluid pressure at the
    // contact interface.  However, since it could be zero,
    // use an average measure of the contact traction instead.
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
    if (bporo) { feLog("    P gap       : %15le", pnorm);
        if (m_atol > 0) feLog("%15le\n", m_atol);
        else feLog("       ***\n");
    }
    for (int isol=0; isol<nsol; ++isol) {
        feLog("    C[%d] gap   : %15le", m_sid[isol], cnorm[isol]);
        if (m_atol > 0) feLog("%15le\n", m_atol);
        else feLog("       ***\n");
    }
    
    feLog("    maximum gap  : %15le", maxgap);
    if (m_gtol > 0) feLog("%15le\n", m_gtol); else feLog("       ***\n");
    if (bporo) {
        feLog("    maximum pgap : %15le", maxpg);
        if (m_ptol > 0) feLog("%15le\n", m_ptol); else feLog("       ***\n");
    }
    for (int isol=0; isol<nsol; ++isol) {
        feLog("    maximum cgap[%d] : %15le", m_sid[isol], maxcg[isol]);
        if (m_ctol > 0) feLog("%15le\n", m_ctol); else feLog("       ***\n");
    }
    
    ProjectSurface(m_ss, m_ms, true);
    if (m_btwo_pass) ProjectSurface(m_ms, m_ss, true);
    
    m_bfreeze = true;
    
    return bconv;
}

//-----------------------------------------------------------------------------
void FESlidingInterfaceMP::Serialize(DumpStream &ar)
{
    // serialize contact data
    FEContactInterface::Serialize(ar);
    
    // serialize contact surface data
    m_ms.Serialize(ar);
    m_ss.Serialize(ar);
    
    
    ar & m_ambp;
    ar & m_ambc;
    ar & m_sid;
    ar & m_ssl;
    ar & m_msl;
    ar & m_sz;
    
    // serialize element pointers
    SerializeElementPointers(m_ss, m_ms, ar);
    SerializeElementPointers(m_ms, m_ss, ar);
}

//-----------------------------------------------------------------------------

void FESlidingInterfaceMP::MarkAmbient()
{	
	int i, j, id, np;
	
    // get number of DOFS
    DOFS& fedofs = GetFEModel()->GetDOFS();
    int MAX_CDOFS = fedofs.GetVariableSize("concentration");
    
	// Mark all nodes as free-draining.  This needs to be done for ALL
	// contact interfaces prior to executing Update(), where nodes that are
	// in contact are subsequently marked as non free-draining.  This ensures
	// that for surfaces involved in more than one contact interface, nodes
	// that have been marked as non free-draining are not reset to 
	// free-draining.
	for (np=0; np<2; ++np)
	{
		FESlidingSurfaceMP& s = (np == 0? m_ss : m_ms);
		
		if (s.m_bporo) {
			// first, mark all nodes as free-draining (= neg. ID)
			// this is done by setting the dof's equation number
			// to a negative number
			for (i=0; i<s.Nodes(); ++i) 
			{
				id = s.Node(i).m_ID[m_dofP];
				if (id >= 0) 
				{
					FENode& node = s.Node(i);
					// mark node as free-draining
					node.m_ID[m_dofP] = -id-2;
				}
			}
		}
		if (s.m_bsolu) {
			// first, mark all nodes as free-draining (= neg. ID)
			// this is done by setting the dof's equation number
			// to a negative number
			for (i=0; i<s.Nodes(); ++i) 
			{
				for (j=0; j<MAX_CDOFS; ++j) {
					id = s.Node(i).m_ID[m_dofC+j];
					if (id >= 0) 
					{
						FENode& node = s.Node(i);
						// mark node as free-draining
						node.m_ID[m_dofC+j] = -id-2;
					}
				}
			}
		}
	}
}

//-----------------------------------------------------------------------------

void FESlidingInterfaceMP::SetAmbient()
{	
	int i, j, np;
	
    // get number of DOFS
    DOFS& fedofs = GetFEModel()->GetDOFS();
    int MAX_CDOFS = fedofs.GetVariableSize("concentration");
    
	// Set the pressure to zero for the free-draining nodes
	for (np=0; np<2; ++np)
	{
		FESlidingSurfaceMP& s = (np == 0? m_ss : m_ms);
		
		if (s.m_bporo) {
			// loop over all nodes
			for (i=0; i<s.Nodes(); ++i) 
			{
				if (s.Node(i).m_ID[m_dofP] < -1)
				{
					FENode& node = s.Node(i);
					// set the fluid pressure to ambient condition
					node.set(m_dofP, m_ambp);
				}
			}
		}
		if (s.m_bsolu) {
			// loop over all nodes
			for (i=0; i<s.Nodes(); ++i) 
			{
				for (j=0; j<MAX_CDOFS; ++j) {
					if (s.Node(i).m_ID[m_dofC+j] < -1)
					{
						FENode& node = s.Node(i);
						// set the fluid pressure to ambient condition
						node.set(m_dofC + j, m_ambc[j]);
					}
				}
			}
		}
	}
}

//-----------------------------------------------------------------------------
FESoluteData* FESlidingInterfaceMP::FindSoluteData(int nid)
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
