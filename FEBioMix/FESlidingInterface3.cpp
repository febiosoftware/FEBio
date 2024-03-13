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
#include "FESlidingInterface3.h"
#include "FEBiphasic.h"
#include "FEBiphasicSolute.h"
#include "FECore/FEModel.h"
#include "FECore/log.h"
#include "FECore/DOFS.h"
#include "FECore/FENormalProjection.h"
#include <FECore/FELinearSystem.h>
#include <FECore/FEAnalysis.h>

//-----------------------------------------------------------------------------
// Define sliding interface parameters
BEGIN_FECORE_CLASS(FESlidingInterface3, FEContactInterface)
	ADD_PARAMETER(m_atol     , "tolerance"            );
    ADD_PARAMETER(m_gtol     , "gaptol"               )->setUnits(UNIT_LENGTH);;
	ADD_PARAMETER(m_ptol     , "ptol"                 );
	ADD_PARAMETER(m_ctol     , "ctol"                 );
	ADD_PARAMETER(m_epsn     , "penalty"              );
	ADD_PARAMETER(m_bautopen , "auto_penalty"         );
    ADD_PARAMETER(m_bupdtpen , "update_penalty"     );
	ADD_PARAMETER(m_btwo_pass, "two_pass"             );
	ADD_PARAMETER(m_knmult   , "knmult"               );
	ADD_PARAMETER(m_stol     , "search_tol"           );
	ADD_PARAMETER(m_epsp     , "pressure_penalty"     );
	ADD_PARAMETER(m_epsc     , "concentration_penalty");
	ADD_PARAMETER(m_bsymm    , "symmetric_stiffness"  );
    ADD_PARAMETER(m_srad     , "search_radius"        )->setUnits(UNIT_LENGTH);;
	ADD_PARAMETER(m_nsegup   , "seg_up"               );
	ADD_PARAMETER(m_naugmin  , "minaug"               );
	ADD_PARAMETER(m_naugmax  , "maxaug"               );
    ADD_PARAMETER(m_breloc   , "node_reloc"         );
    ADD_PARAMETER(m_bsmaug   , "smooth_aug"         );
	ADD_PARAMETER(m_ambp     , "ambient_pressure"     );
	ADD_PARAMETER(m_ambc     , "ambient_concentration");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FESlidingSurface3::Data::Data()
{
	m_nu = vec3d(0,0,0);
	m_rs = vec2d(0,0);
	m_Lmd = 0.0;
	m_Lmc = 0.0;
	m_Lmp = 0.0;
	m_epsn = 1.0;
	m_epsc = 1.0;
	m_epsp = 1.0;
	m_pg = 0.0;
	m_cg = 0.0;
}

//-----------------------------------------------------------------------------
void FESlidingSurface3::Data::Serialize(DumpStream& ar)
{
	FEBiphasicContactPoint::Serialize(ar);
	ar & m_Lmd;
	ar & m_Lmc;
	ar & m_epsn;
	ar & m_epsp;
	ar & m_epsc;
	ar & m_cg;
	ar & m_nu;
	ar & m_rs;
}

//-----------------------------------------------------------------------------
// FESlidingSurface3
//-----------------------------------------------------------------------------

FESlidingSurface3::FESlidingSurface3(FEModel* pfem) : FEBiphasicContactSurface(pfem)
{ 
	m_bporo = m_bsolu = false;
	m_dofC = -1;
}

//-----------------------------------------------------------------------------
FESlidingSurface3::~FESlidingSurface3()
{

}

//-----------------------------------------------------------------------------
//! create material point data
FEMaterialPoint* FESlidingSurface3::CreateMaterialPoint()
{
	return new FESlidingSurface3::Data;
}

//-----------------------------------------------------------------------------
void FESlidingSurface3::UnpackLM(FEElement& el, vector<int>& lm)
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
bool FESlidingSurface3::Init()
{
	// initialize surface data first
	if (FEBiphasicContactSurface::Init() == false) return false;
	
	// store concentration index
	DOFS& dofs = GetFEModel()->GetDOFS();
	m_dofC = dofs.GetDOF("concentration", 0);

	// allocate data structures
	m_nn.assign(Nodes(), vec3d(0,0,0));
    m_pn.assign(Nodes(), 0);
	
	// determine biphasic and biphasic-solute status
	m_poro.resize(Elements(),false);
	m_solu.resize(Elements(),-1);
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
			FEBiphasic* pb = dynamic_cast<FEBiphasic*> (pm);
			FEBiphasicSolute* pbs = dynamic_cast<FEBiphasicSolute*> (pm);
			if (pb || pbs) {
				m_poro[i] = true;
				m_bporo = true;
			}
			if (pbs) {
				m_solu[i] = pbs->GetSolute()->GetSoluteID() - 1;
				m_bsolu = true;
			}
		}
	}

	return true;
}

//-----------------------------------------------------------------------------
//! Evaluate the nodal contact pressures by averaging values from surrounding
//! faces.  This function ensures that nodal contact pressures are always
//! positive, so that they can be used to detect free-draining status.

void FESlidingSurface3::EvaluateNodalContactPressures()
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
//! This function calculates the node normal. Due to the piecewise continuity
//! of the surface elements this normal is not uniquely defined so in order to
//! obtain a unique normal the normal is averaged for each node over all the 
//! element normals at the node

void FESlidingSurface3::UpdateNodeNormals()
{
	vec3d y[FEElement::MAX_NODES];
	
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
	int N = Nodes();
	for (int i=0; i<N; ++i) m_nn[i].unit();
}

//-----------------------------------------------------------------------------
vec3d FESlidingSurface3::GetContactForce()
{
    return m_Ft;
}

//-----------------------------------------------------------------------------
double FESlidingSurface3::GetContactArea()
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
vec3d FESlidingSurface3::GetFluidForce()
{
	int n, i;
	const int MN = FEElement::MAX_NODES;
	double pn[MN];
    
    // get parent contact interface to extract ambient fluid pressure
    FESlidingInterface3* si3 = dynamic_cast<FESlidingInterface3*>(GetContactInterface());
    double ambp = si3->m_ambp;
	
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
void FESlidingSurface3::Serialize(DumpStream& ar)
{
	FEBiphasicContactSurface::Serialize(ar);
	ar & m_dofC;
	ar & m_dofP;
	ar & m_bporo;
	ar & m_bsolu;
	ar & m_nn;
	ar & m_poro;
	ar & m_solu;
	ar & m_pn;
	ar & m_Ft;
}

//-----------------------------------------------------------------------------
void FESlidingSurface3::GetContactPressure(int nface, double& pg)
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
void FESlidingSurface3::GetContactTraction(int nface, vec3d& pt)
{
    FESurfaceElement& el = Element(nface);
    int ni = el.GaussPoints();
    pt = vec3d(0,0,0);
	for (int k = 0; k < ni; ++k)
	{
		Data& data = static_cast<Data&>(*el.GetMaterialPoint(k));
		pt -= data.m_nu*data.m_Ln;
	}
    pt /= ni;
}

//-----------------------------------------------------------------------------
void FESlidingSurface3::GetNodalContactPressure(int nface, double* pn)
{
	FESurfaceElement& el = Element(nface);
    for (int k=0; k<el.Nodes(); ++k)
        pn[k] = m_pn[el.m_lnode[k]];
}

//-----------------------------------------------------------------------------
void FESlidingSurface3::GetNodalContactTraction(int nface, vec3d* tn)
{
	FESurfaceElement& el = Element(nface);
    for (int k=0; k<el.Nodes(); ++k)
        tn[k] = m_nn[el.m_lnode[k]]*(-m_pn[el.m_lnode[k]]);
}

//-----------------------------------------------------------------------------
// FESlidingInterface3
//-----------------------------------------------------------------------------

FESlidingInterface3::FESlidingInterface3(FEModel* pfem) : FEContactInterface(pfem), m_ss(pfem), m_ms(pfem)
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
	m_ambc = 0;
	m_nsegup = 0;
	m_bautopen = false;
    m_bupdtpen = false;
    m_breloc = false;
    m_bsmaug = false;
	
	m_naugmin = 0;
	m_naugmax = 10;

	m_dofP = pfem->GetDOFIndex("p");
	m_dofC = pfem->GetDOFIndex("concentration", 0);
	
	m_ss.SetSibling(&m_ms);
	m_ms.SetSibling(&m_ss);
    m_ss.SetContactInterface(this);
    m_ms.SetContactInterface(this);
}

//-----------------------------------------------------------------------------

FESlidingInterface3::~FESlidingInterface3()
{
}

//-----------------------------------------------------------------------------
bool FESlidingInterface3::Init()
{
	m_Rgas = GetFEModel()->GetGlobalConstant("R");
	m_Tabs = GetFEModel()->GetGlobalConstant("T");

	// initialize surface data
	if (m_ss.Init() == false) return false;
	if (m_ms.Init() == false) return false;
	
	return true;
}

//-----------------------------------------------------------------------------
//! build the matrix profile for use in the stiffness matrix
void FESlidingInterface3::BuildMatrixProfile(FEGlobalMatrix& K)
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
	const int dof_C = fem.GetDOFIndex("concentration", 0);

	vector<int> lm(8*FEElement::MAX_NODES*2);
					
	int npass = (m_btwo_pass?2:1);
	for (int np=0; np<npass; ++np)
	{
		FESlidingSurface3& ss = (np == 0? m_ss : m_ms);
		FESlidingSurface3& ms = (np == 0? m_ms : m_ss);
		// get the mesh
		int sid, mid;
						
		int k, l;
		for (int j=0; j<ss.Elements(); ++j)
		{
			FESurfaceElement& se = ss.Element(j);
			sid = ss.m_solu[j];
			int nint = se.GaussPoints();
			int* sn = &se.m_node[0];
			for (k=0; k<nint; ++k)
			{
				FESlidingSurface3::Data& data = static_cast<FESlidingSurface3::Data&>(*se.GetMaterialPoint(k));

				FESurfaceElement* pe = data.m_pme;
				if (pe != 0)
				{
					FESurfaceElement& me = *pe;
					mid = ms.m_solu[pe->m_lid];
					int* mn = &me.m_node[0];
									
					assign(lm, -1);
									
					int nseln = se.Nodes();
					int nmeln = me.Nodes();
									
					for (l=0; l<nseln; ++l)
					{
						vector<int>& id = mesh.Node(sn[l]).m_ID;
						lm[8*l  ] = id[dof_X];
						lm[8*l+1] = id[dof_Y];
						lm[8*l+2] = id[dof_Z];
						lm[8*l+3] = id[dof_P];
						lm[8*l+4] = id[dof_RU];
						lm[8*l+5] = id[dof_RV];
						lm[8*l+6] = id[dof_RW];
						lm[8*l+7] = id[dof_C + sid];
					}
									
					for (l=0; l<nmeln; ++l)
					{
						vector<int>& id = mesh.Node(mn[l]).m_ID;
						lm[8*(l+nseln)  ] = id[dof_X];
						lm[8*(l+nseln)+1] = id[dof_Y];
						lm[8*(l+nseln)+2] = id[dof_Z];
						lm[8*(l+nseln)+3] = id[dof_P];
						lm[8*(l+nseln)+4] = id[dof_RU];
						lm[8*(l+nseln)+5] = id[dof_RV];
						lm[8*(l+nseln)+6] = id[dof_RW];
						lm[8*(l+nseln)+7] = id[dof_C + mid];
					}
									
					K.build_add(lm);
				}
			}
		}
	}

}

//-----------------------------------------------------------------------------
void FESlidingInterface3::UpdateAutoPenalty()
{
    // calculate the penalty
    if (m_bautopen)
    {
        CalcAutoPenalty(m_ss);
        CalcAutoPenalty(m_ms);
        if (m_ss.m_bporo) CalcAutoPressurePenalty(m_ss);
        if (m_ss.m_bsolu) CalcAutoConcentrationPenalty(m_ss);
        if (m_ms.m_bporo) CalcAutoPressurePenalty(m_ms);
        if (m_ms.m_bsolu) CalcAutoConcentrationPenalty(m_ms);
    }
}

//-----------------------------------------------------------------------------
void FESlidingInterface3::Activate()
{
	// don't forget to call the base class
	FEContactInterface::Activate();

    UpdateAutoPenalty();
    
	// update sliding interface data
	Update();
}

//-----------------------------------------------------------------------------
void FESlidingInterface3::CalcAutoPenalty(FESlidingSurface3& s)
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
			FESlidingSurface3::Data& pt = static_cast<FESlidingSurface3::Data&>(*el.GetMaterialPoint(j));
			pt.m_epsn = eps;
		}
	}
}

//-----------------------------------------------------------------------------
void FESlidingInterface3::CalcAutoPressurePenalty(FESlidingSurface3& s)
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
			FESlidingSurface3::Data& pt = static_cast<FESlidingSurface3::Data&>(*el.GetMaterialPoint(j));
			pt.m_epsp = eps;
		}
	}
}

//-----------------------------------------------------------------------------
double FESlidingInterface3::AutoPressurePenalty(FESurfaceElement& el, FESlidingSurface3& s)
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
	double eps = 0.0;
	FEBiphasic* bp = dynamic_cast<FEBiphasic*> (pm);
	FEBiphasicSolute* bps = dynamic_cast<FEBiphasicSolute*> (pm);
	if (bp)
	{
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
		bp->Permeability(K, mp);
			
		eps = n.x*(K[0][0]*n.x+K[0][1]*n.y+K[0][2]*n.z)
		+n.y*(K[1][0]*n.x+K[1][1]*n.y+K[1][2]*n.z)
		+n.z*(K[2][0]*n.x+K[2][1]*n.y+K[2][2]*n.z);
	}
	else if (bps)
	{
		// get a material point
		FEMaterialPoint& mp = *pe->GetMaterialPoint(0);
		FEElasticMaterialPoint& ept = *(mp.ExtractData<FEElasticMaterialPoint>());
			
		// setup the material point
		ept.m_F = mat3dd(1.0);
		ept.m_J = 1;
		ept.m_s.zero();
			
		// if this is a biphasic-solute element, then get the permeability tensor
		FEBiphasicMaterialPoint& ppt = *(mp.ExtractData<FEBiphasicMaterialPoint>());
		FESolutesMaterialPoint& spt = *(mp.ExtractData<FESolutesMaterialPoint>());
		ppt.m_p = 0;
		ppt.m_w = vec3d(0,0,0);
		spt.m_c[0] = 0;
		spt.m_j[0] = vec3d(0,0,0);
			
		mat3ds K = bps->GetPermeability()->Permeability(mp);
			
		eps = n*(K*n);
	}

	// get the area of the surface element
	double A = s.FaceArea(el);

	// get the volume of the volume element
	double V = m.ElementVolume(*pe);

	return eps*A/V;
}

//-----------------------------------------------------------------------------
void FESlidingInterface3::CalcAutoConcentrationPenalty(FESlidingSurface3& s)
{
	// loop over all surface elements
	for (int i=0; i<s.Elements(); ++i)
	{
		// get the surface element
		FESurfaceElement& el = s.Element(i);
		
		// calculate a modulus
		double eps = AutoConcentrationPenalty(el, s);
		
		// assign to integation points of surface element
		int nint = el.GaussPoints();
		for (int j=0; j<nint; ++j)
		{
			FESlidingSurface3::Data& pt = static_cast<FESlidingSurface3::Data&>(*el.GetMaterialPoint(j));
			pt.m_epsc = eps;
		}
	}
}

//-----------------------------------------------------------------------------
double FESlidingInterface3::AutoConcentrationPenalty(FESurfaceElement& el, FESlidingSurface3& s)
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
		
	// see if this is a biphasic-solute element
	FEBiphasicSolute* pbs = dynamic_cast<FEBiphasicSolute*> (pm);
	if (pbs == 0) return 0.0;

	// get a material point
	FEMaterialPoint& mp = *pe->GetMaterialPoint(0);
	FEElasticMaterialPoint& ept = *(mp.ExtractData<FEElasticMaterialPoint>());
			
	// setup the material point
	ept.m_F = mat3dd(1.0);
	ept.m_J = 1;
	ept.m_s.zero();
			
	// if this is a biphasic-solute element, then get the diffusivity tensor
	FEBiphasicMaterialPoint& ppt = *(mp.ExtractData<FEBiphasicMaterialPoint>());
	FESolutesMaterialPoint& spt = *(mp.ExtractData<FESolutesMaterialPoint>());
	ppt.m_p = 0;
	ppt.m_w = vec3d(0,0,0);
	spt.m_c[0] = 0;
	spt.m_j[0] = vec3d(0,0,0);
			
	mat3ds D = pbs->GetSolute()->m_pDiff->Diffusivity(mp)
	*(pbs->Porosity(mp)*pbs->GetSolute()->m_pSolub->Solubility(mp));
			
	double eps = n*(D*n);
	
	// get the area of the surface element
	double A = s.FaceArea(el);

	// get the volume of the volume element
	double V = m.ElementVolume(*pe);

	return eps*A/V;
}

//-----------------------------------------------------------------------------
void FESlidingInterface3::ProjectSurface(FESlidingSurface3& ss, FESlidingSurface3& ms, bool bupseg, bool bmove)
{
	FEMesh& mesh = GetFEModel()->GetMesh();
	FESurfaceElement* pme;
	vec3d r, nu;
	double rs[2];
	double Ln;
	
	double ps[FEElement::MAX_NODES], p1;
	double cs[FEElement::MAX_NODES], c1;
	
	double R = m_srad*mesh.GetBoundingBox().radius();
	
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
//    #pragma omp parallel for shared(R, bupseg)
	for (int i=0; i<ss.Elements(); ++i)
	{
		FESurfaceElement& el = ss.Element(i);

		bool sporo = ss.m_poro[i];
		int sid = ss.m_solu[i];
		bool ssolu = (sid > -1) ? true : false;
		
		int ne = el.Nodes();
		int nint = el.GaussPoints();
		
		// get the nodal pressures
		if (sporo)
		{
			for (int j=0; j<ne; ++j) ps[j] = mesh.Node(el.m_node[j]).get(m_dofP);
		}
		
		// get the nodal concentrations
		if (ssolu)
		{
			for (int j=0; j<ne; ++j) cs[j] = mesh.Node(el.m_node[j]).get(m_dofC + sid);
		}
		
		for (int j=0; j<nint; ++j)
		{
			FESlidingSurface3::Data& pt = static_cast<FESlidingSurface3::Data&>(*el.GetMaterialPoint(j));

			// calculate the global position of the integration point
			r = ss.Local2Global(el, j);
			
			// get the pressure at the integration point
			if (sporo) p1 = el.eval(ps, j);
			
			// get the concentration at the integration point
			if (ssolu) c1 = el.eval(cs, j);
			
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
				
				pt.m_gap = (g <= R? g : 0);
				
				if ((Ln >= 0) && (g <= R))
				{
					
					// calculate the pressure gap function
					bool mporo = ms.m_poro[pme->m_lid];
					int mid = ms.m_solu[pme->m_lid];
					bool msolu = (mid > -1) ? true : false;
					if (sporo && mporo) {
						double pm[FEElement::MAX_NODES];
						for (int k=0; k<pme->Nodes(); ++k) pm[k] = mesh.Node(pme->m_node[k]).get(m_dofP);
						double p2 = pme->eval(pm, rs[0], rs[1]);
						pt.m_pg = p1 - p2;
					}
					if (ssolu && msolu && (sid == mid)) {
						double cm[FEElement::MAX_NODES];
						for (int k=0; k<pme->Nodes(); ++k) cm[k] = mesh.Node(pme->m_node[k]).get(m_dofC + sid);
						double c2 = pme->eval(cm, rs[0], rs[1]);
						pt.m_cg = c1 - c2;
					}
				}
				else
				{
                    pt.m_Lmd = 0;
					pt.m_gap = 0;
					pt.m_pme = 0;
					if (sporo) {
                        pt.m_Lmp = 0;
						pt.m_pg = 0;
					}
					if (ssolu) {
                        pt.m_Lmc = 0;
						pt.m_cg = 0;
					}
				}
			}
			else
			{
				// the node is not in contact
				pt.m_Lmd = 0;
				pt.m_gap = 0;
				if (sporo) {
					pt.m_Lmp = 0;
					pt.m_pg = 0;
				}
				if (ssolu) {
					pt.m_Lmc = 0;
					pt.m_cg = 0;
				}
			}
		}
	}
}

//-----------------------------------------------------------------------------

void FESlidingInterface3::Update()
{	
	double rs[2];

	FEModel& fem = *GetFEModel();

    // get number of DOFS
    DOFS& fedofs = GetFEModel()->GetDOFS();
    int MAX_CDOFS = fedofs.GetVariableSize("concentration");
    
	double R = m_srad*fem.GetMesh().GetBoundingBox().radius();
	
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
	//	Logfile& log = GetLogfile();
	//	log.printf("seg_up iteration # %d\n", niter+1);
	
	// project the surfaces onto each other
	// this will update the gap functions as well
    static bool bfirst = true;
    ProjectSurface(m_ss, m_ms, bupseg, (m_breloc && bfirst));
	if (m_btwo_pass || m_ss.m_bporo) ProjectSurface(m_ms, m_ss, bupseg);
    bfirst = false;
	
	// Update the net contact pressures
	UpdateContactPressures();

    // update node normals
    m_ss.UpdateNodeNormals();
    m_ms.UpdateNodeNormals();
    
	// set poro flag
	bool bporo = (m_ss.m_bporo || m_ms.m_bporo);
	
	// only continue if we are doing a poro-elastic simulation
	if (bporo == false) return;
	
	// Now that the nodes have been projected, we need to figure out
	// if we need to modify the constraints on the pressure and concentration dofs.
	// If the nodes are not in contact, they must be match ambient
	// conditions. Since all nodes have been previously marked to be
	// under ambient conditions in MarkAmbient(), we just need to reverse
	// this setting here, for nodes that are in contact.
	
	// Next, we loop over each surface, visiting the nodes
	// and finding out if that node is in contact or not
	int npass = (m_btwo_pass?2:1);
	for (int np=0; np<npass; ++np)
	{
		FESlidingSurface3& ss = (np == 0? m_ss : m_ms);
		FESlidingSurface3& ms = (np == 0? m_ms : m_ss);
		
        // loop over all the nodes of the primary surface
        for (int n=0; n<ss.Nodes(); ++n) {
            if (ss.m_pn[n] > 0) {
                FENode& node = ss.Node(n);
                int id = node.m_ID[m_dofP];
                if (id < -1)
                    // mark node as non-free-draining (= pos ID)
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
            // initialize projection data
            FENormalProjection np(ss);
            np.SetTolerance(m_stol);
            np.SetSearchRadius(m_srad);
            np.Init();
            
			for (int n = 0; n<ms.Nodes(); ++n)
			{
				// get the node
				FENode& node = ms.Node(n);
				
				// project it onto the primary surface
				FESurfaceElement* pse = np.Project(node.m_rt, ms.m_nn[n], rs);
				
				if (pse)
				{
					// we found an element, so let's see if it's even remotely close to contact
					// find the global location of the intersection point
					vec3d q = ss.Local2Global(*pse, rs[0], rs[1]);
					
					// calculate the gap function
					double g = ms.m_nn[n]*(node.m_rt - q);
					
					if (fabs(g) <= R)
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
							for (int j = 0; j<MAX_CDOFS; ++j) {
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
void FESlidingInterface3::LoadVector(FEGlobalVector& R, const FETimeInfo& tp)
{
	vector<int> sLM, mLM, LM, en;
	vector<double> fe;
	const int MN = FEElement::MAX_NODES;
	double detJ[MN], w[MN], *Hs, Hm[MN];
	double N[10*MN];

	FEModel& fem = *GetFEModel();

	double dt = fem.GetTime().timeIncrement;
	
    double psf = GetPenaltyScaleFactor();
    
    m_ss.m_Ft = vec3d(0, 0, 0);
    m_ms.m_Ft = vec3d(0, 0, 0);
    
	// loop over the nr of passes
	int npass = (m_btwo_pass?2:1);
	for (int np=0; np<npass; ++np)
	{
		// get primary and secondary surface
		FESlidingSurface3& ss = (np == 0? m_ss : m_ms);
		FESlidingSurface3& ms = (np == 0? m_ms : m_ss);
		
		// loop over all primary surface elements
		for (int i = 0; i<ss.Elements(); ++i)
		{
			// get the surface element
			FESurfaceElement& se = ss.Element(i);

			bool sporo = ss.m_poro[i];
			int sid = ss.m_solu[i];
			bool ssolu = (sid > -1) ? true : false;
			
			// get the nr of nodes and integration points
			int nseln = se.Nodes();
			int nint = se.GaussPoints();
			
			// copy the LM vector; we'll need it later
			ss.UnpackLM(se, sLM);
			
			// we calculate all the metrics we need before we
			// calculate the nodal forces
			for (int j = 0; j<nint; ++j)
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
			for (int j = 0; j<nint; ++j)
			{
				FESlidingSurface3::Data& pt = static_cast<FESlidingSurface3::Data&>(*se.GetMaterialPoint(j));
				// get the secondary surface element
				FESurfaceElement* pme = pt.m_pme;
				if (pme)
				{
					// get the secondary surface element
					FESurfaceElement& me = *pme;

					bool mporo = ms.m_poro[pme->m_lid];
					int mid = ms.m_solu[pme->m_lid];
					bool msolu = (mid > -1) ? true : false;
					
					// get the nr of secondary surface element nodes
					int nmeln = me.Nodes();
					
					// copy LM vector
					ms.UnpackLM(me, mLM);
					
					// calculate degrees of freedom
					int ndof = 3*(nseln + nmeln);
					
					// build the LM vector
					LM.resize(ndof);
					for (int k = 0; k<nseln; ++k)
					{
						LM[3*k  ] = sLM[3*k  ];
						LM[3*k+1] = sLM[3*k+1];
						LM[3*k+2] = sLM[3*k+2];
					}
					
					for (int k = 0; k<nmeln; ++k)
					{
						LM[3*(k+nseln)  ] = mLM[3*k  ];
						LM[3*(k+nseln)+1] = mLM[3*k+1];
						LM[3*(k+nseln)+2] = mLM[3*k+2];
					}
					
					// build the en vector
					en.resize(nseln+nmeln);
					for (int k = 0; k<nseln; ++k) en[k] = se.m_node[k];
					for (int k = 0; k<nmeln; ++k) en[k + nseln] = me.m_node[k];
					
					// get element shape functions
					Hs = se.H(j);
					
					// get secondary surface element shape functions
					double r = pt.m_rs[0];
					double s = pt.m_rs[1];
					me.shape_fnc(Hm, r, s);
					
					// get normal vector
					vec3d nu = pt.m_nu;
					
					// gap function
					double g = pt.m_gap;
					
					// lagrange multiplier
					double Lm = pt.m_Lmd;
					
					// penalty 
					double eps = m_epsn*pt.m_epsn*psf;
					
					// contact traction
					double tn = Lm + eps*g;
					tn = MBRACKET(tn);
					
					// calculate the force vector
					fe.resize(ndof);
					zero(fe);
					
					for (int k = 0; k<nseln; ++k)
					{
						N[3*k  ] = -Hs[k]*nu.x;
						N[3*k+1] = -Hs[k]*nu.y;
						N[3*k+2] = -Hs[k]*nu.z;
					}
					
					for (int k = 0; k<nmeln; ++k)
					{
						N[3*(k+nseln)  ] = Hm[k]*nu.x;
						N[3*(k+nseln)+1] = Hm[k]*nu.y;
						N[3*(k+nseln)+2] = Hm[k]*nu.z;
					}
					
					for (int k = 0; k<ndof; ++k) fe[k] += tn*N[k] * detJ[j] * w[j];
					
                    for (int k=0; k<nseln; ++k)
                    {
                        ss.m_Ft += vec3d(fe[k*3], fe[k*3+1], fe[k*3+2]);
                    }
                    
                    for (int k = 0; k<nmeln; ++k)
                    {
                        ms.m_Ft += vec3d(fe[(k + nseln) * 3], fe[(k + nseln) * 3 + 1], fe[(k + nseln) * 3 + 2]);
                    }
                    
					// assemble the global residual
					R.Assemble(en, LM, fe);
					
					// do the biphasic stuff
					if (tn > 0)
					{
						if (sporo && mporo)
						{
							// calculate nr of pressure dofs
							int ndof = nseln + nmeln;
							
							// calculate the flow rate
							double epsp = m_epsp*pt.m_epsp*psf;
							
							double wn = pt.m_Lmp + epsp*pt.m_pg;
							
							// fill the LM
							LM.resize(ndof);
							for (int k = 0; k<nseln; ++k) LM[k] = sLM[3 * nseln + k];
							for (int k = 0; k<nmeln; ++k) LM[k + nseln] = mLM[3 * nmeln + k];
							
							// fill the force array
							fe.resize(ndof);
							zero(fe);
							for (int k = 0; k<nseln; ++k) N[k] = Hs[k];
							for (int k = 0; k<nmeln; ++k) N[k + nseln] = -Hm[k];
							
							for (int k = 0; k<ndof; ++k) fe[k] += dt*wn*N[k] * detJ[j] * w[j];
							
							// assemble residual
							R.Assemble(en, LM, fe);
						}
						if (ssolu && msolu && (sid == mid))
						{
							// calculate nr of concentration dofs
							int ndof = nseln + nmeln;
							
							// calculate the flow rate
							double epsc = m_epsc*pt.m_epsc*psf;
							
							double jn = pt.m_Lmc + epsc*pt.m_cg;
							
							// fill the LM
							LM.resize(ndof);
							for (int k=0; k<nseln; ++k) LM[k        ] = sLM[(4+sid)*nseln+k];
							for (int k=0; k<nmeln; ++k) LM[k + nseln] = mLM[(4+mid)*nmeln+k];
							
							// fill the force array
							fe.resize(ndof);
							zero(fe);
							for (int k=0; k<nseln; ++k) N[k      ] =  Hs[k];
							for (int k=0; k<nmeln; ++k) N[k+nseln] = -Hm[k];
							
							for (int k = 0; k<ndof; ++k) fe[k] += dt*jn*N[k] * detJ[j] * w[j];
							
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
void FESlidingInterface3::StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp)
{
	int i, j, k, l;
	vector<int> sLM, mLM, LM, en;
	const int MN = FEElement::MAX_NODES;
	double detJ[MN], w[MN], *Hs, Hm[MN];
	double pt[MN], dpr[MN], dps[MN];
	double ct[MN], dcr[MN], dcs[MN];
	double N[10*MN];
	FEElementMatrix ke;

	FEModel& fem = *GetFEModel();

    double psf = GetPenaltyScaleFactor();
    
	// see how many reformations we've had to do so far
	int nref = LS.GetSolver()->m_nref;
	
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
	
	// do single- or two-pass
	int npass = (m_btwo_pass?2:1);
	for (int np=0; np < npass; ++np)
	{
		// get the primary and secondary surface
		FESlidingSurface3& ss = (np == 0? m_ss : m_ms);
		FESlidingSurface3& ms = (np == 0? m_ms : m_ss);
		
		// loop over all primary surface elements
		for (i=0; i<ss.Elements(); ++i)
		{
			// get the next element
			FESurfaceElement& se = ss.Element(i);

			bool sporo = ss.m_poro[i];
			int sid = ss.m_solu[i];
			bool ssolu = (sid > -1) ? true : false;
			
			// get nr of nodes and integration points
			int nseln = se.Nodes();
			int nint = se.GaussPoints();

			double pn[FEElement::MAX_NODES], cn[FEElement::MAX_NODES];
			for (j=0; j<nseln; ++j)
			{
				pn[j] = ss.GetMesh()->Node(se.m_node[j]).get(m_dofP);
				cn[j] = ss.GetMesh()->Node(se.m_node[j]).get(m_dofC + sid);
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
				
				// pressure
				if (sporo)
				{
					pt[j] = se.eval(pn, j);
					dpr[j] = se.eval_deriv1(pn, j);
					dps[j] = se.eval_deriv2(pn, j);
				}
				// concentration
				if (ssolu)
				{
					ct[j] = se.eval(cn, j);
					dcr[j] = se.eval_deriv1(cn, j);
					dcs[j] = se.eval_deriv2(cn, j);
				}
			}
			
			// loop over all integration points
			for (j=0; j<nint; ++j)
			{
				FESlidingSurface3::Data& pt = static_cast<FESlidingSurface3::Data&>(*se.GetMaterialPoint(j));

				// get the secondary surface element
				FESurfaceElement* pme = pt.m_pme;
				if (pme)
				{
					FESurfaceElement& me = *pme;

					bool mporo = ms.m_poro[pme->m_lid];
					int mid = ms.m_solu[pme->m_lid];
					bool msolu = (mid > -1) ? true : false;
					
					// get the nr of nodes
					int nmeln = me.Nodes();

					// nodal data
					double pm[FEElement::MAX_NODES], cm[FEElement::MAX_NODES];
					for (k=0; k<nmeln; ++k)
					{
						pm[k] = ms.GetMesh()->Node(me.m_node[k]).get(m_dofP);
						cm[k] = ms.GetMesh()->Node(me.m_node[k]).get(m_dofC + mid);
					}
					
					// copy the LM vector
					ms.UnpackLM(me, mLM);
					
					int ndpn;	// number of dofs per node
					int ndof;	// number of dofs in stiffness matrix
					
					if (ssolu && msolu && (sid == mid)) {
						// calculate dofs for biphasic-solute contact
						ndpn = 5;
						ndof = ndpn*(nseln+nmeln);
						
						// build the LM vector
						LM.resize(ndof);
						
						for (k=0; k<nseln; ++k)
						{
							LM[ndpn*k  ] = sLM[3*k  ];			// x-dof
							LM[ndpn*k+1] = sLM[3*k+1];			// y-dof
							LM[ndpn*k+2] = sLM[3*k+2];			// z-dof
							LM[ndpn*k+3] = sLM[3*nseln+k];		// p-dof
							LM[ndpn*k+4] = sLM[(4+sid)*nseln+k];		// c-dof
						}
						for (k=0; k<nmeln; ++k)
						{
							LM[ndpn*(k+nseln)  ] = mLM[3*k  ];			// x-dof
							LM[ndpn*(k+nseln)+1] = mLM[3*k+1];			// y-dof
							LM[ndpn*(k+nseln)+2] = mLM[3*k+2];			// z-dof
							LM[ndpn*(k+nseln)+3] = mLM[3*nmeln+k];		// p-dof
							LM[ndpn*(k+nseln)+4] = mLM[(4+mid)*nmeln+k];		// c-dof
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
					
					// gap function
					double g = pt.m_gap;
					
					// lagrange multiplier
					double Lm = pt.m_Lmd;
					
					// penalty 
					double eps = m_epsn*pt.m_epsn*psf;
					
					// contact traction
					double tn = Lm + eps*g;
					tn = MBRACKET(tn);
					
					//	double dtn = m_eps*HEAVYSIDE(Lm + eps*g);
					double dtn = (tn > 0.? eps :0.);
					
					// create the stiffness matrix
					ke.resize(ndof, ndof); ke.zero();
					
					// --- S O L I D - S O L I D   C O N T A C T ---
					
					// a. NxN-term
					//------------------------------------
					
					// calculate the N-vector
					for (k=0; k<nseln; ++k)
					{
						N[ndpn*k  ] = Hs[k]*nu.x;
						N[ndpn*k+1] = Hs[k]*nu.y;
						N[ndpn*k+2] = Hs[k]*nu.z;
					}
					
					for (k=0; k<nmeln; ++k)
					{
						N[ndpn*(k+nseln)  ] = -Hm[k]*nu.x;
						N[ndpn*(k+nseln)+1] = -Hm[k]*nu.y;
						N[ndpn*(k+nseln)+2] = -Hm[k]*nu.z;
					}
					
					if (ndpn == 5) {
						for (k=0; k<nseln; ++k)
						{
							N[ndpn*k+3] = 0;
							N[ndpn*k+4] = 0;
						}
						for (k=0; k<nmeln; ++k)
						{
							N[ndpn*(k+nseln)+3] = 0;
							N[ndpn*(k+nseln)+4] = 0;
						}
					}
					else if (ndpn == 4) {
						for (k=0; k<nseln; ++k)
							N[ndpn*k+3] = 0;
						for (k=0; k<nmeln; ++k)
							N[ndpn*(k+nseln)+3] = 0;
					}
					
					for (k=0; k<ndof; ++k)
						for (l=0; l<ndof; ++l) ke[k][l] += dtn*N[k]*N[l]*detJ[j]*w[j];
					
					// b. A-term
					//-------------------------------------
					
					for (k=0; k<nseln; ++k) N[k      ] =  Hs[k];
					for (k=0; k<nmeln; ++k) N[k+nseln] = -Hm[k];
					
					double* Gr = se.Gr(j);
					double* Gs = se.Gs(j);
					vec3d gs[2];
					ss.CoBaseVectors(se, j, gs);
					
					mat3d S1, S2;
					S1.skew(gs[0]);
					S2.skew(gs[1]);
					mat3d As[FEElement::MAX_NODES];
					for (l=0; l<nseln; ++l)
						As[l] = S2*Gr[l] - S1*Gs[l];
					
					if (!m_bsymm)
					{	// non-symmetric
						for (l=0; l<nseln; ++l)
						{
							for (k=0; k<nseln+nmeln; ++k)
							{
								ke[k*ndpn  ][l*ndpn  ] -= knmult*tn*w[j]*N[k]*As[l][0][0];
								ke[k*ndpn  ][l*ndpn+1] -= knmult*tn*w[j]*N[k]*As[l][0][1];
								ke[k*ndpn  ][l*ndpn+2] -= knmult*tn*w[j]*N[k]*As[l][0][2];
								
								ke[k*ndpn+1][l*ndpn  ] -= knmult*tn*w[j]*N[k]*As[l][1][0];
								ke[k*ndpn+1][l*ndpn+1] -= knmult*tn*w[j]*N[k]*As[l][1][1];
								ke[k*ndpn+1][l*ndpn+2] -= knmult*tn*w[j]*N[k]*As[l][1][2];
								
								ke[k*ndpn+2][l*ndpn  ] -= knmult*tn*w[j]*N[k]*As[l][2][0];
								ke[k*ndpn+2][l*ndpn+1] -= knmult*tn*w[j]*N[k]*As[l][2][1];
								ke[k*ndpn+2][l*ndpn+2] -= knmult*tn*w[j]*N[k]*As[l][2][2];
							}
						}
					} 
					else 
					{	// symmetric
						for (l=0; l<nseln; ++l)
						{
							for (k=0; k<nseln+nmeln; ++k)
							{
								ke[k*ndpn  ][l*ndpn  ] -= 0.5*knmult*tn*w[j]*N[k]*As[l][0][0];
								ke[k*ndpn  ][l*ndpn+1] -= 0.5*knmult*tn*w[j]*N[k]*As[l][0][1];
								ke[k*ndpn  ][l*ndpn+2] -= 0.5*knmult*tn*w[j]*N[k]*As[l][0][2];
								
								ke[k*ndpn+1][l*ndpn  ] -= 0.5*knmult*tn*w[j]*N[k]*As[l][1][0];
								ke[k*ndpn+1][l*ndpn+1] -= 0.5*knmult*tn*w[j]*N[k]*As[l][1][1];
								ke[k*ndpn+1][l*ndpn+2] -= 0.5*knmult*tn*w[j]*N[k]*As[l][1][2];
								
								ke[k*ndpn+2][l*ndpn  ] -= 0.5*knmult*tn*w[j]*N[k]*As[l][2][0];
								ke[k*ndpn+2][l*ndpn+1] -= 0.5*knmult*tn*w[j]*N[k]*As[l][2][1];
								ke[k*ndpn+2][l*ndpn+2] -= 0.5*knmult*tn*w[j]*N[k]*As[l][2][2];
								
								ke[l*ndpn  ][k*ndpn  ] -= 0.5*knmult*tn*w[j]*N[k]*As[l][0][0];
								ke[l*ndpn+1][k*ndpn  ] -= 0.5*knmult*tn*w[j]*N[k]*As[l][0][1];
								ke[l*ndpn+2][k*ndpn  ] -= 0.5*knmult*tn*w[j]*N[k]*As[l][0][2];
								
								ke[l*ndpn  ][k*ndpn+1] -= 0.5*knmult*tn*w[j]*N[k]*As[l][1][0];
								ke[l*ndpn+1][k*ndpn+1] -= 0.5*knmult*tn*w[j]*N[k]*As[l][1][1];
								ke[l*ndpn+2][k*ndpn+1] -= 0.5*knmult*tn*w[j]*N[k]*As[l][1][2];
								
								ke[l*ndpn  ][k*ndpn+2] -= 0.5*knmult*tn*w[j]*N[k]*As[l][2][0];
								ke[l*ndpn+1][k*ndpn+2] -= 0.5*knmult*tn*w[j]*N[k]*As[l][2][1];
								ke[l*ndpn+2][k*ndpn+2] -= 0.5*knmult*tn*w[j]*N[k]*As[l][2][2];
							}
						}
					}
					
					// c. M-term
					//---------------------------------------
					
					vec3d Gm[2];
					ms.ContraBaseVectors(me, r, s, Gm);
					
					// evaluate secondary surface normal
					vec3d mnu = Gm[0] ^ Gm[1];
					mnu.unit();
					
					double Hmr[FEElement::MAX_NODES], Hms[FEElement::MAX_NODES];
					me.shape_deriv(Hmr, Hms, r, s);
					vec3d mm[FEElement::MAX_NODES];
					for (k=0; k<nmeln; ++k) 
						mm[k] = Gm[0]*Hmr[k] + Gm[1]*Hms[k];
					
					if (!m_bsymm)
					{	// non-symmetric
						for (k=0; k<nmeln; ++k) 
						{
							for (l=0; l<nseln+nmeln; ++l)
							{
								ke[(k+nseln)*ndpn  ][l*ndpn  ] += tn*knmult*detJ[j]*w[j]*mnu.x*mm[k].x*N[l];
								ke[(k+nseln)*ndpn  ][l*ndpn+1] += tn*knmult*detJ[j]*w[j]*mnu.x*mm[k].y*N[l];
								ke[(k+nseln)*ndpn  ][l*ndpn+2] += tn*knmult*detJ[j]*w[j]*mnu.x*mm[k].z*N[l];
								
								ke[(k+nseln)*ndpn+1][l*ndpn  ] += tn*knmult*detJ[j]*w[j]*mnu.y*mm[k].x*N[l];
								ke[(k+nseln)*ndpn+1][l*ndpn+1] += tn*knmult*detJ[j]*w[j]*mnu.y*mm[k].y*N[l];
								ke[(k+nseln)*ndpn+1][l*ndpn+2] += tn*knmult*detJ[j]*w[j]*mnu.y*mm[k].z*N[l];
								
								ke[(k+nseln)*ndpn+2][l*ndpn  ] += tn*knmult*detJ[j]*w[j]*mnu.z*mm[k].x*N[l];
								ke[(k+nseln)*ndpn+2][l*ndpn+1] += tn*knmult*detJ[j]*w[j]*mnu.z*mm[k].y*N[l];
								ke[(k+nseln)*ndpn+2][l*ndpn+2] += tn*knmult*detJ[j]*w[j]*mnu.z*mm[k].z*N[l];
							}
						}
					}
					else
					{	// symmetric
						for (k=0; k<nmeln; ++k) 
						{
							for (l=0; l<nseln+nmeln; ++l)
							{
								ke[(k+nseln)*ndpn  ][l*ndpn  ] += 0.5*knmult*tn*detJ[j]*w[j]*mnu.x*mm[k].x*N[l];
								ke[(k+nseln)*ndpn  ][l*ndpn+1] += 0.5*knmult*tn*detJ[j]*w[j]*mnu.x*mm[k].y*N[l];
								ke[(k+nseln)*ndpn  ][l*ndpn+2] += 0.5*knmult*tn*detJ[j]*w[j]*mnu.x*mm[k].z*N[l];
								
								ke[(k+nseln)*ndpn+1][l*ndpn  ] += 0.5*knmult*tn*detJ[j]*w[j]*mnu.y*mm[k].x*N[l];
								ke[(k+nseln)*ndpn+1][l*ndpn+1] += 0.5*knmult*tn*detJ[j]*w[j]*mnu.y*mm[k].y*N[l];
								ke[(k+nseln)*ndpn+1][l*ndpn+2] += 0.5*knmult*tn*detJ[j]*w[j]*mnu.y*mm[k].z*N[l];
								
								ke[(k+nseln)*ndpn+2][l*ndpn  ] += 0.5*knmult*tn*detJ[j]*w[j]*mnu.z*mm[k].x*N[l];
								ke[(k+nseln)*ndpn+2][l*ndpn+1] += 0.5*knmult*tn*detJ[j]*w[j]*mnu.z*mm[k].y*N[l];
								ke[(k+nseln)*ndpn+2][l*ndpn+2] += 0.5*knmult*tn*detJ[j]*w[j]*mnu.z*mm[k].z*N[l];
								
								ke[l*ndpn  ][(k+nseln)*ndpn  ] += 0.5*knmult*tn*detJ[j]*w[j]*mnu.x*mm[k].x*N[l];
								ke[l*ndpn+1][(k+nseln)*ndpn  ] += 0.5*knmult*tn*detJ[j]*w[j]*mnu.y*mm[k].x*N[l];
								ke[l*ndpn+2][(k+nseln)*ndpn  ] += 0.5*knmult*tn*detJ[j]*w[j]*mnu.z*mm[k].x*N[l];
								
								ke[l*ndpn  ][(k+nseln)*ndpn+1] += 0.5*knmult*tn*detJ[j]*w[j]*mnu.x*mm[k].y*N[l];
								ke[l*ndpn+1][(k+nseln)*ndpn+1] += 0.5*knmult*tn*detJ[j]*w[j]*mnu.y*mm[k].y*N[l];
								ke[l*ndpn+2][(k+nseln)*ndpn+1] += 0.5*knmult*tn*detJ[j]*w[j]*mnu.z*mm[k].y*N[l];
								
								ke[l*ndpn  ][(k+nseln)*ndpn+2] += 0.5*knmult*tn*detJ[j]*w[j]*mnu.x*mm[k].z*N[l];
								ke[l*ndpn+1][(k+nseln)*ndpn+2] += 0.5*knmult*tn*detJ[j]*w[j]*mnu.y*mm[k].z*N[l];
								ke[l*ndpn+2][(k+nseln)*ndpn+2] += 0.5*knmult*tn*detJ[j]*w[j]*mnu.z*mm[k].z*N[l];
							}
						}
					}
					
					// --- B I P H A S I C - S O L U T E  S T I F F N E S S ---
					if (ssolu && msolu && (sid == mid))
					{
						double dt = fem.GetTime().timeIncrement;
						
						double epsp = (tn > 0) ? m_epsp*pt.m_epsp*psf : 0.;
						double epsc = (tn > 0) ? m_epsc*pt.m_epsc*psf : 0.;
						
						// --- S O L I D - P R E S S U R E / S O L U T E   C O N T A C T ---
						
						if (!m_bsymm)
						{
							
							// a. q-term
							//-------------------------------------
							
							double dpmr, dpms;
							dpmr = me.eval_deriv1(pm, r, s);
							dpms = me.eval_deriv2(pm, r, s);
							
							double dcmr, dcms;
							dcmr = me.eval_deriv1(cm, r, s);
							dcms = me.eval_deriv2(cm, r, s);
							
							for (k=0; k<nseln+nmeln; ++k)
								for (l=0; l<nseln+nmeln; ++l)
								{
									ke[ndpn*k + 3][ndpn*l  ] += dt*w[j]*detJ[j]*epsp*N[k]*N[l]*(dpmr*Gm[0].x + dpms*Gm[1].x);
									ke[ndpn*k + 3][ndpn*l+1] += dt*w[j]*detJ[j]*epsp*N[k]*N[l]*(dpmr*Gm[0].y + dpms*Gm[1].y);
									ke[ndpn*k + 3][ndpn*l+2] += dt*w[j]*detJ[j]*epsp*N[k]*N[l]*(dpmr*Gm[0].z + dpms*Gm[1].z);

									ke[ndpn*k + 4][ndpn*l  ] += dt*w[j]*detJ[j]*epsc*N[k]*N[l]*(dcmr*Gm[0].x + dcms*Gm[1].x);
									ke[ndpn*k + 4][ndpn*l+1] += dt*w[j]*detJ[j]*epsc*N[k]*N[l]*(dcmr*Gm[0].y + dcms*Gm[1].y);
									ke[ndpn*k + 4][ndpn*l+2] += dt*w[j]*detJ[j]*epsc*N[k]*N[l]*(dcmr*Gm[0].z + dcms*Gm[1].z);
								}
							
							double wn = pt.m_Lmp + epsp*pt.m_pg;
							double jn = pt.m_Lmc + epsc*pt.m_cg;
							
							// b. A-term
							//-------------------------------------
							
							for (l=0; l<nseln; ++l)
								for (k=0; k<nseln+nmeln; ++k)
								{
									ke[ndpn*k + 3][ndpn*l  ] -= dt*w[j]*wn*N[k]*(As[l][0][0]*nu.x + As[l][0][1]*nu.y + As[l][0][2]*nu.z);
									ke[ndpn*k + 3][ndpn*l+1] -= dt*w[j]*wn*N[k]*(As[l][1][0]*nu.x + As[l][1][1]*nu.y + As[l][1][2]*nu.z);
									ke[ndpn*k + 3][ndpn*l+2] -= dt*w[j]*wn*N[k]*(As[l][2][0]*nu.x + As[l][2][1]*nu.y + As[l][2][2]*nu.z);

									ke[ndpn*k + 4][ndpn*l  ] -= dt*w[j]*jn*N[k]*(As[l][0][0]*nu.x + As[l][0][1]*nu.y + As[l][0][2]*nu.z);
									ke[ndpn*k + 4][ndpn*l+1] -= dt*w[j]*jn*N[k]*(As[l][1][0]*nu.x + As[l][1][1]*nu.y + As[l][1][2]*nu.z);
									ke[ndpn*k + 4][ndpn*l+2] -= dt*w[j]*jn*N[k]*(As[l][2][0]*nu.x + As[l][2][1]*nu.y + As[l][2][2]*nu.z);
								}
							
							// c. m-term
							//---------------------------------------
							
							for (k=0; k<nmeln; ++k)
								for (l=0; l<nseln+nmeln; ++l)
								{
									ke[ndpn*(k+nseln) + 3][ndpn*l  ] += dt*w[j]*detJ[j]*wn*N[l]*mm[k].x;
									ke[ndpn*(k+nseln) + 3][ndpn*l+1] += dt*w[j]*detJ[j]*wn*N[l]*mm[k].y;
									ke[ndpn*(k+nseln) + 3][ndpn*l+2] += dt*w[j]*detJ[j]*wn*N[l]*mm[k].z;

									ke[ndpn*(k+nseln) + 4][ndpn*l  ] += dt*w[j]*detJ[j]*jn*N[l]*mm[k].x;
									ke[ndpn*(k+nseln) + 4][ndpn*l+1] += dt*w[j]*detJ[j]*jn*N[l]*mm[k].y;
									ke[ndpn*(k+nseln) + 4][ndpn*l+2] += dt*w[j]*detJ[j]*jn*N[l]*mm[k].z;
								}
						}
						
						
						// --- P R E S S U R E - P R E S S U R E   C O N T A C T ---
						
						// calculate the N-vector
						for (k=0; k<nseln; ++k)
						{
							N[ndpn*k+3] = Hs[k];
							N[ndpn*k+4] = Hs[k];
						}
						
						for (k=0; k<nmeln; ++k)
						{
							N[ndpn*(k+nseln)+3] = -Hm[k];
							N[ndpn*(k+nseln)+4] = -Hm[k];
						}
						
						for (k=3; k<ndof; k+=ndpn)
							for (l=3; l<ndof; l+=ndpn) ke[k][l] -= dt*epsp*w[j]*detJ[j]*N[k]*N[l];

						// --- C O N C E N T R A T I O N - C O N C E N T R A T I O N   C O N T A C T ---
						
						for (k=4; k<ndof; k+=ndpn)
							for (l=4; l<ndof; l+=ndpn) ke[k][l] -= dt*epsc*w[j]*detJ[j]*N[k]*N[l];
						
					}
					
					// --- B I P H A S I C   S T I F F N E S S ---
					else if (sporo && mporo)
					{
						// the variable dt is either the timestep or one
						// depending on whether we are using the symmetric
						// poro version or not.
						double dt = fem.GetTime().timeIncrement;
						
						double epsp = (tn > 0) ? m_epsp*pt.m_epsp*psf : 0.;
						
						// --- S O L I D - P R E S S U R E   C O N T A C T ---
						
						if (!m_bsymm)
						{
							
							// a. q-term
							//-------------------------------------
							
							double dpmr, dpms;
							dpmr = me.eval_deriv1(pm, r, s);
							dpms = me.eval_deriv2(pm, r, s);
							
							for (k=0; k<nseln+nmeln; ++k)
								for (l=0; l<nseln+nmeln; ++l)
								{
									ke[ndpn*k + 3][ndpn*l  ] += dt*w[j]*detJ[j]*epsp*N[k]*N[l]*(dpmr*Gm[0].x + dpms*Gm[1].x);
									ke[ndpn*k + 3][ndpn*l+1] += dt*w[j]*detJ[j]*epsp*N[k]*N[l]*(dpmr*Gm[0].y + dpms*Gm[1].y);
									ke[ndpn*k + 3][ndpn*l+2] += dt*w[j]*detJ[j]*epsp*N[k]*N[l]*(dpmr*Gm[0].z + dpms*Gm[1].z);
								}
							
							double wn = pt.m_Lmp + epsp*pt.m_pg;
							
							// b. A-term
							//-------------------------------------
							
							for (l=0; l<nseln; ++l)
								for (k=0; k<nseln+nmeln; ++k)
								{
									ke[ndpn*k + 3][ndpn*l  ] -= dt*w[j]*wn*N[k]*(As[l][0][0]*nu.x + As[l][0][1]*nu.y + As[l][0][2]*nu.z);
									ke[ndpn*k + 3][ndpn*l+1] -= dt*w[j]*wn*N[k]*(As[l][1][0]*nu.x + As[l][1][1]*nu.y + As[l][1][2]*nu.z);
									ke[ndpn*k + 3][ndpn*l+2] -= dt*w[j]*wn*N[k]*(As[l][2][0]*nu.x + As[l][2][1]*nu.y + As[l][2][2]*nu.z);
								}
							
							// c. m-term
							//---------------------------------------
							
							for (k=0; k<nmeln; ++k)
								for (l=0; l<nseln+nmeln; ++l)
								{
									ke[ndpn*(k+nseln) + 3][ndpn*l  ] += dt*w[j]*detJ[j]*wn*N[l]*mm[k].x;
									ke[ndpn*(k+nseln) + 3][ndpn*l+1] += dt*w[j]*detJ[j]*wn*N[l]*mm[k].y;
									ke[ndpn*(k+nseln) + 3][ndpn*l+2] += dt*w[j]*detJ[j]*wn*N[l]*mm[k].z;
								}
						}
						
						
						// --- P R E S S U R E - P R E S S U R E   C O N T A C T ---
						
						// calculate the N-vector
						for (k=0; k<nseln; ++k)
						{
							N[ndpn*k+3] = Hs[k];
						}
						
						for (k=0; k<nmeln; ++k)
						{
							N[ndpn*(k+nseln)+3] = -Hm[k];
						}
						
						for (k=3; k<ndof; k+=ndpn)
							for (l=3; l<ndof; l+=ndpn) ke[k][l] -= dt*epsp*w[j]*detJ[j]*N[k]*N[l];
						
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
void FESlidingInterface3::UpdateContactPressures()
{
	int npass = (m_btwo_pass?2:1);
	for (int np=0; np<npass; ++np)
	{
		FESlidingSurface3& ss = (np == 0? m_ss : m_ms);
		FESlidingSurface3& ms = (np == 0? m_ms : m_ss);
		
		// loop over all elements of the primary surface
		for (int n=0; n<ss.Elements(); ++n)
		{
			FESurfaceElement& el = ss.Element(n);
			int nint = el.GaussPoints();
			
			// get the normal tractions at the integration points
			double gap, eps;
			for (int i=0; i<nint; ++i) 
			{
				FESlidingSurface3::Data& pt = static_cast<FESlidingSurface3::Data&>(*el.GetMaterialPoint(i));

				gap = pt.m_gap;
				eps = m_epsn*pt.m_epsn;
				pt.m_Ln = MBRACKET(pt.m_Lmd + eps*gap);
				FESurfaceElement* pme = pt.m_pme;
				if (m_btwo_pass && pme)
				{
					int mint = pme->GaussPoints();
					double ti[FEElement::MAX_NODES];
					for (int j=0; j<mint; ++j) {
						FESlidingSurface3::Data& md = static_cast<FESlidingSurface3::Data&>(*pme->GetMaterialPoint(j));

						gap = md.m_gap;
						eps = m_epsn*md.m_epsn;
						ti[j] = MBRACKET(md.m_Lmd + m_epsn*md.m_epsn*md.m_gap);
					}
					// project the data to the nodes
					double tn[FEElement::MAX_NODES];
					pme->FEElement::project_to_nodes(ti, tn);
					// now evaluate the traction at the intersection point
					double Ln = pme->eval(tn, pt.m_rs[0], pt.m_rs[1]);
					pt.m_Ln += MBRACKET(Ln);
				}
			}
		}
        ss.EvaluateNodalContactPressures();
	}
}

//-----------------------------------------------------------------------------
bool FESlidingInterface3::Augment(int naug, const FETimeInfo& tp)
{
	// make sure we need to augment
	if (m_laugon != 1) return true;

	double Ln, Lp, Lc;
	bool bconv = true;
	
	bool bporo = (m_ss.m_bporo && m_ms.m_bporo);
	bool bsolu = (m_ss.m_bsolu && m_ms.m_bsolu);

	int NS = m_ss.Elements();
	int NM = m_ms.Elements();
	
	// --- c a l c u l a t e   i n i t i a l   n o r m s ---
	// a. normal component
	double normL0 = 0, normP = 0, normDP = 0, normC = 0, normDC = 0;
	for (int i=0; i<NS; ++i)
	{
		FESurfaceElement& el = m_ss.Element(i);
		for (int j=0; j<el.GaussPoints(); ++j)
		{
			FESlidingSurface3::Data& ds = static_cast<FESlidingSurface3::Data&>(*el.GetMaterialPoint(j));
			normL0 += ds.m_Lmd*ds.m_Lmd;
		}
	}
	for (int i=0; i<NM; ++i)
	{
		FESurfaceElement& el = m_ms.Element(i);
		for (int j=0; j<el.GaussPoints(); ++j)
		{
			FESlidingSurface3::Data& dm = static_cast<FESlidingSurface3::Data&>(*el.GetMaterialPoint(j));
			normL0 += dm.m_Lmd*dm.m_Lmd;
		}
	}
	
	// b. gap component
	// (is calculated during update)
	double maxgap = 0, maxpg = 0, maxcg = 0;
	
	// update Lagrange multipliers
	double normL1 = 0, eps, epsp, epsc;
    for (int i=0; i<m_ss.Elements(); ++i) {
        FESurfaceElement& el = m_ss.Element(i);
        vec3d tn[FEElement::MAX_INTPOINTS];
        if (m_bsmaug) m_ss.GetGPSurfaceTraction(i, tn);
        for (int j=0; j<el.GaussPoints(); ++j) {
			FESlidingSurface3::Data& data = static_cast<FESlidingSurface3::Data&>(*el.GetMaterialPoint(j));
			// update Lagrange multipliers on primary surface
            if (m_bsmaug) {
                // replace this multiplier with a smoother version
                Ln = -(tn[j]*data.m_nu);
                data.m_Lmd = MBRACKET(Ln);
                if (m_btwo_pass) data.m_Lmd /= 2;
            }
            else {
                eps = m_epsn*data.m_epsn;
                Ln = data.m_Lmd + eps*data.m_gap;
                data.m_Lmd = MBRACKET(Ln);
            }
            
            normL1 += data.m_Lmd*data.m_Lmd;
            
            if (m_ss.m_bporo) {
                Lp = 0;
                if (Ln > 0) {
                    epsp = m_epsp*data.m_epsp;
                    Lp = data.m_Lmp + epsp*data.m_pg;
                    maxpg = max(maxpg,fabs(data.m_pg));
                    normDP += data.m_pg*data.m_pg;
                }
                data.m_Lmp = Lp;
            }
            
            if (m_ss.m_bsolu) {
                Lc = 0;
                if (Ln > 0) {
                    epsc = m_epsc*data.m_epsc;
                    Lc = data.m_Lmc + epsc*data.m_cg;
                    maxcg = max(maxcg,fabs(data.m_cg));
                    normDC += data.m_cg*data.m_cg;
                }
                data.m_Lmc = Lc;
            }
            
            if (Ln > 0) maxgap = max(maxgap,fabs(data.m_gap));
        }
    }
    
    for (int i=0; i<m_ms.Elements(); ++i) {
        FESurfaceElement& el = m_ms.Element(i);
        vec3d tn[FEElement::MAX_INTPOINTS];
        if (m_bsmaug) m_ms.GetGPSurfaceTraction(i, tn);
        for (int j=0; j<el.GaussPoints(); ++j) {
			FESlidingSurface3::Data& data = static_cast<FESlidingSurface3::Data&>(*el.GetMaterialPoint(j));
			// update Lagrange multipliers on secondary surface
            if (m_bsmaug) {
                // replace this multiplier with a smoother version
                Ln = -(tn[j]*data.m_nu);
                data.m_Lmd = MBRACKET(Ln);
                if (m_btwo_pass) data.m_Lmd /= 2;
            }
            else {
                eps = m_epsn*data.m_epsn;
                Ln = data.m_Lmd + eps*data.m_gap;
                data.m_Lmd = MBRACKET(Ln);
            }
            
            normL1 += data.m_Lmd*data.m_Lmd;
            
            if (m_ms.m_bporo) {
                Lp = 0;
                if (Ln > 0) {
                    epsp = m_epsp*data.m_epsp;
                    Lp = data.m_Lmp + epsp*data.m_pg;
                    maxpg = max(maxpg,fabs(data.m_pg));
                    normDP += data.m_pg*data.m_pg;
                }
                data.m_Lmp = Lp;
            }
            
            if (m_ms.m_bsolu) {
                Lc = 0;
                if (Ln > 0) {
                    epsc = m_epsc*data.m_epsc;
                    Lc = data.m_Lmc + epsc*data.m_cg;
                    maxcg = max(maxcg,fabs(data.m_cg));
                    normDC += data.m_cg*data.m_cg;
                }
                data.m_Lmc = Lc;
            }
            
            if (Ln > 0) maxgap = max(maxgap,fabs(data.m_gap));
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
	double cnorm = (normC != 0 ? (normDC/normC) : normDC); 
	
	// check convergence
	if ((m_gtol > 0) && (maxgap > m_gtol)) bconv = false;
	if ((m_ptol > 0) && (bporo && maxpg > m_ptol)) bconv = false;
	if ((m_ctol > 0) && (bsolu && maxcg > m_ctol)) bconv = false;
	
	if ((m_atol > 0) && (lnorm > m_atol)) bconv = false;
	if ((m_atol > 0) && (pnorm > m_atol)) bconv = false;
	if ((m_atol > 0) && (cnorm > m_atol)) bconv = false;

	if (naug < m_naugmin ) bconv = false;
	if (naug >= m_naugmax) bconv = true;

	feLog(" sliding interface # %d\n", GetID());
	feLog("                        CURRENT        REQUIRED\n");
	feLog("    D multiplier : %15le", lnorm); if (m_atol > 0) feLog("%15le\n", m_atol); else feLog("       ***\n");
	if (bporo) { feLog("    P gap        : %15le", pnorm); if (m_atol > 0) feLog("%15le\n", m_atol); else feLog("       ***\n"); }
	if (bsolu) { feLog("    C gap        : %15le", cnorm); if (m_atol > 0) feLog("%15le\n", m_atol); else feLog("       ***\n"); }
	
	feLog("    maximum gap  : %15le", maxgap);
	if (m_gtol > 0) feLog("%15le\n", m_gtol); else feLog("       ***\n");
	if (bporo) {
		feLog("    maximum pgap : %15le", maxpg);
		if (m_ptol > 0) feLog("%15le\n", m_ptol); else feLog("       ***\n");
	}
	if (bsolu) {
		feLog("    maximum cgap : %15le", maxcg);
		if (m_ctol > 0) feLog("%15le\n", m_ctol); else feLog("       ***\n");
	}
	
	return bconv;
}

//-----------------------------------------------------------------------------
void FESlidingInterface3::Serialize(DumpStream &ar)
{
	// store contact data
	FEContactInterface::Serialize(ar);

	// store contact surface data
	m_ms.Serialize(ar);
	m_ss.Serialize(ar);

	// serialize element pointers
	SerializeElementPointers(m_ss, m_ms, ar);
	SerializeElementPointers(m_ms, m_ss, ar);
}

//-----------------------------------------------------------------------------

void FESlidingInterface3::MarkAmbient()
{	
	// Mark all nodes as free-draining.  This needs to be done for ALL
	// contact interfaces prior to executing Update(), where nodes that are
	// in contact are subsequently marked as non free-draining.  This ensures
	// that for surfaces involved in more than one contact interface, nodes
	// that have been marked as non free-draining are not reset to 
	// free-draining.

    // get number of DOFS
    DOFS& fedofs = GetFEModel()->GetDOFS();
    int MAX_CDOFS = fedofs.GetVariableSize("concentration");
    
	for (int np=0; np<2; ++np)
	{
		FESlidingSurface3& s = (np == 0? m_ss : m_ms);
		
		if (s.m_bporo) {
			// first, mark all nodes as free-draining (= neg. ID)
			// this is done by setting the dof's equation number
			// to a negative number
			for (int i=0; i<s.Nodes(); ++i) 
			{
				int id = s.Node(i).m_ID[m_dofP];
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
			for (int i=0; i<s.Nodes(); ++i) 
			{
				for (int j=0; j<MAX_CDOFS; ++j) {
					int id = s.Node(i).m_ID[m_dofC+j];
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

void FESlidingInterface3::SetAmbient()
{	
    // get number of DOFS
    DOFS& fedofs = GetFEModel()->GetDOFS();
    int MAX_CDOFS = fedofs.GetVariableSize("concentration");
    
	// Set the pressure to zero for the free-draining nodes
	for (int np=0; np<2; ++np)
	{
		FESlidingSurface3& s = (np == 0? m_ss : m_ms);
		
		if (s.m_bporo) {
			// loop over all nodes
			for (int i=0; i<s.Nodes(); ++i) 
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
			for (int i=0; i<s.Nodes(); ++i) 
			{
				for (int j=0; j<MAX_CDOFS; ++j)
				{
					if (s.Node(i).m_ID[m_dofC+j] < -1)
					{
						FENode& node = s.Node(i);
						// set the fluid pressure to ambient condition
						node.set(m_dofC + j, m_ambc);
					}
				}
			}
		}
	}
}
