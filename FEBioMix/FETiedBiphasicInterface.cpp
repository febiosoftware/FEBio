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
#include "FETiedBiphasicInterface.h"
#include "FEBiphasic.h"
#include "FECore/FEModel.h"
#include "FECore/FEAnalysis.h"
#include "FECore/FENormalProjection.h"
#include <FECore/FELinearSystem.h>
#include "FECore/log.h"

//-----------------------------------------------------------------------------
// Define sliding interface parameters
BEGIN_FECORE_CLASS(FETiedBiphasicInterface, FEContactInterface)
	ADD_PARAMETER(m_atol     , "tolerance"          );
	ADD_PARAMETER(m_gtol     , "gaptol"             );
	ADD_PARAMETER(m_ptol     , "ptol"               );
	ADD_PARAMETER(m_epsn     , "penalty"            );
	ADD_PARAMETER(m_bautopen , "auto_penalty"       );
    ADD_PARAMETER(m_bupdtpen , "update_penalty"     );
	ADD_PARAMETER(m_btwo_pass, "two_pass"           );
	ADD_PARAMETER(m_knmult   , "knmult"             );
	ADD_PARAMETER(m_stol     , "search_tol"         );
	ADD_PARAMETER(m_epsp     , "pressure_penalty"   );
	ADD_PARAMETER(m_bsymm    , "symmetric_stiffness");
	ADD_PARAMETER(m_srad     , "search_radius"      );
	ADD_PARAMETER(m_naugmin  , "minaug"             );
	ADD_PARAMETER(m_naugmax  , "maxaug"             );
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FETiedBiphasicSurface::Data::Data()
{
    m_Gap = vec3d(0,0,0);
    m_dg = vec3d(0,0,0);
    m_nu = vec3d(0,0,0);
    m_rs = vec2d(0,0);
    m_Lmd = vec3d(0,0,0);
    m_tr = vec3d(0,0,0);
    m_Lmp  = 0.0;
    m_epsn = 1.0;
    m_epsp = 1.0;
    m_pg   = 0.0;
}

//-----------------------------------------------------------------------------
void FETiedBiphasicSurface::Data::Serialize(DumpStream& ar)
{
	FEBiphasicContactPoint::Serialize(ar);
	ar & m_Gap;
	ar & m_dg;
	ar & m_nu;
	ar & m_rs;
	ar & m_Lmd;
	ar & m_tr;
	ar & m_epsn;
	ar & m_epsp;
}

//-----------------------------------------------------------------------------
// FETiedBiphasicSurface
//-----------------------------------------------------------------------------

FETiedBiphasicSurface::FETiedBiphasicSurface(FEModel* pfem) : FEBiphasicContactSurface(pfem)
{ 
	m_bporo = false;
}

//-----------------------------------------------------------------------------
//! create material point data
FEMaterialPoint* FETiedBiphasicSurface::CreateMaterialPoint()
{
	return new FETiedBiphasicSurface::Data;
}

//-----------------------------------------------------------------------------
bool FETiedBiphasicSurface::Init()
{
	// initialize surface data first
	if (FEBiphasicContactSurface::Init() == false) return false;
	
    // allocate node normals
    m_nn.assign(Nodes(), vec3d(0,0,0));
    
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
//! This function calculates the node normal. Due to the piecewise continuity
//! of the surface elements this normal is not uniquely defined so in order to
//! obtain a unique normal the normal is averaged for each node over all the 
//! element normals at the node

void FETiedBiphasicSurface::UpdateNodeNormals()
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
void FETiedBiphasicSurface::Serialize(DumpStream& ar)
{
	FEBiphasicContactSurface::Serialize(ar);
	ar & m_bporo;
	ar & m_poro;
	ar & m_nn;
}

//-----------------------------------------------------------------------------
void FETiedBiphasicSurface::GetVectorGap(int nface, vec3d& pg)
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
void FETiedBiphasicSurface::GetContactTraction(int nface, vec3d& pt)
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
// FETiedBiphasicInterface
//-----------------------------------------------------------------------------

FETiedBiphasicInterface::FETiedBiphasicInterface(FEModel* pfem) : FEContactInterface(pfem), m_ss(pfem), m_ms(pfem)
{
	static int count = 1;
	SetID(count++);
	
	// initial values
	m_knmult = 1;
	m_atol = 0.1;
	m_epsn = 1;
	m_epsp = 1;
	m_btwo_pass = false;
	m_stol = 0.01;
	m_bsymm = true;
	m_srad = 1.0;
	m_gtol = -1;	// we use augmentation tolerance by default
	m_ptol = -1;	// we use augmentation tolerance by default
	m_bautopen = false;
    m_bupdtpen = false;
	
	m_naugmin = 0;
	m_naugmax = 10;

	m_dofP = pfem->GetDOFIndex("p");
	
	// set parents
	m_ss.SetContactInterface(this);
	m_ms.SetContactInterface(this);

	m_ss.SetSibling(&m_ms);
	m_ms.SetSibling(&m_ss);
}

//-----------------------------------------------------------------------------

FETiedBiphasicInterface::~FETiedBiphasicInterface()
{
}

//-----------------------------------------------------------------------------
bool FETiedBiphasicInterface::Init()
{
	// initialize surface data
	if (m_ss.Init() == false) return false;
	if (m_ms.Init() == false) return false;

	return true;
}

//-----------------------------------------------------------------------------
//! build the matrix profile for use in the stiffness matrix
void FETiedBiphasicInterface::BuildMatrixProfile(FEGlobalMatrix& K)
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
		FETiedBiphasicSurface& ss = (np == 0? m_ss : m_ms);
		FETiedBiphasicSurface& ms = (np == 0? m_ms : m_ss);
						
		int ni = 0, k, l;
		for (int j=0; j<ss.Elements(); ++j)
		{
			FESurfaceElement& se = ss.Element(j);
			int nint = se.GaussPoints();
			int* sn = &se.m_node[0];
			for (k=0; k<nint; ++k, ++ni)
			{
				FETiedBiphasicSurface::Data& pt = static_cast<FETiedBiphasicSurface::Data&>(*se.GetMaterialPoint(k));
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
void FETiedBiphasicInterface::UpdateAutoPenalty()
{
    // calculate the penalty
    if (m_bautopen)
    {
        CalcAutoPenalty(m_ss);
        if (m_ss.m_bporo) CalcAutoPressurePenalty(m_ss);
        if (m_btwo_pass) {
            CalcAutoPenalty(m_ms);
            if (m_ms.m_bporo) CalcAutoPressurePenalty(m_ms);
        }
    }
}

//-----------------------------------------------------------------------------
void FETiedBiphasicInterface::Activate()
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
void FETiedBiphasicInterface::CalcAutoPenalty(FETiedBiphasicSurface& s)
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
			FETiedBiphasicSurface::Data& pt = static_cast<FETiedBiphasicSurface::Data&>(*el.GetMaterialPoint(j));
			pt.m_epsn = eps;
        }
	}
}

//-----------------------------------------------------------------------------
void FETiedBiphasicInterface::CalcAutoPressurePenalty(FETiedBiphasicSurface& s)
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
			FETiedBiphasicSurface::Data& pt = static_cast<FETiedBiphasicSurface::Data&>(*el.GetMaterialPoint(j));
			pt.m_epsp = eps;
        }
	}
}

//-----------------------------------------------------------------------------
double FETiedBiphasicInterface::AutoPressurePenalty(FESurfaceElement& el, FETiedBiphasicSurface& s)
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
// Perform initial projection between tied surfaces in reference configuration
void FETiedBiphasicInterface::InitialProjection(FETiedBiphasicSurface& ss, FETiedBiphasicSurface& ms)
{
	FEMesh& mesh = GetFEModel()->GetMesh();
	double R = m_srad*mesh.GetBoundingBox().radius();
	
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
			
			FETiedBiphasicSurface::Data& pt = static_cast<FETiedBiphasicSurface::Data&>(*el.GetMaterialPoint(j));
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
void FETiedBiphasicInterface::ProjectSurface(FETiedBiphasicSurface& ss, FETiedBiphasicSurface& ms)
{
	FEMesh& mesh = GetFEModel()->GetMesh();
	FESurfaceElement* pme;
	vec3d r;
	
	double ps[FEElement::MAX_NODES], p1;
	
	// loop over all integration points
	for (int i=0; i<ss.Elements(); ++i)
	{
		FESurfaceElement& el = ss.Element(i);
		bool sporo = ss.m_poro[i];
		
		int ne = el.Nodes();
		int nint = el.GaussPoints();
		
		// get the nodal pressures
		if (sporo)
		{
			for (int j=0; j<ne; ++j) ps[j] = mesh.Node(el.m_node[j]).get(m_dofP);
		}
		
		for (int j=0; j<nint; ++j)
		{
			FETiedBiphasicSurface::Data& pt = static_cast<FETiedBiphasicSurface::Data&>(*el.GetMaterialPoint(j));

			// calculate the global position of the integration point
			r = ss.Local2Global(el, j);
			
			// get the pressure at the integration point
			if (sporo) p1 = el.eval(ps, j);
			
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
			}
			else
			{
				// the node is not tied
				pt.m_dg = vec3d(0,0,0);
				if (sporo) pt.m_pg = 0;
			}
		}
	}
}

//-----------------------------------------------------------------------------

void FETiedBiphasicInterface::Update()
{	
	// project the surfaces onto each other
	// this will update the gap functions as well
	ProjectSurface(m_ss, m_ms);
	if (m_btwo_pass) ProjectSurface(m_ms, m_ss);
	
}

//-----------------------------------------------------------------------------
void FETiedBiphasicInterface::LoadVector(FEGlobalVector& R, const FETimeInfo& tp)
{
	int i, j, k;
	vector<int> sLM, mLM, LM, en;
	vector<double> fe;
	const int MN = FEElement::MAX_NODES;
	double detJ[MN], w[MN], *Hs, Hm[MN];
	double N[8*MN];

	// get time step
	// if we're using the symmetric formulation
	// we need to multiply with the timestep
    double dt = GetFEModel()->GetTime().timeIncrement;
	
    // Update auto-penalty if requested
    if (m_bupdtpen && GetFEModel()->GetCurrentStep()->GetFESolver()->m_niter == 0) UpdateAutoPenalty();
    
	// loop over the nr of passes
	int npass = (m_btwo_pass?2:1);
	for (int np=0; np<npass; ++np)
	{
		// get primary and secondary surface
		FETiedBiphasicSurface& ss = (np == 0? m_ss : m_ms);
		FETiedBiphasicSurface& ms = (np == 0? m_ms : m_ss);
		
		// loop over all primary surface elements
		for (i=0; i<ss.Elements(); ++i)
		{
			// get the surface element
			FESurfaceElement& se = ss.Element(i);
			
			bool sporo = ss.m_poro[i];
			
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
			}
			
			// loop over all integration points
			// note that we are integrating over the current surface
			for (j=0; j<nint; ++j)
			{
				FETiedBiphasicSurface::Data& pt = static_cast<FETiedBiphasicSurface::Data&>(*se.GetMaterialPoint(j));

				// get the secondary surface element
				FESurfaceElement* pme = pt.m_pme;
				if (pme)
				{
					// get the secondary surface element
					FESurfaceElement& me = *pme;
					
					bool mporo = ms.m_poro[pme->m_lid];
					
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
					
					// get normal vector
					vec3d nu = pt.m_nu;
					
					// gap function
					vec3d dg = pt.m_dg;
					
					// lagrange multiplier
					vec3d Lm = pt.m_Lmd;
					
					// penalty 
					double eps = m_epsn*pt.m_epsn;
					
					// contact traction
					vec3d t = Lm + dg*eps;
                    pt.m_tr = t;
					
					// calculate the force vector
					fe.resize(ndof);
					zero(fe);
					
					for (k=0; k<nseln; ++k)
					{
						N[3*k  ] = Hs[k]*t.x;
						N[3*k+1] = Hs[k]*t.y;
						N[3*k+2] = Hs[k]*t.z;
					}
					
					for (k=0; k<nmeln; ++k)
					{
						N[3*(k+nseln)  ] = -Hm[k]*t.x;
						N[3*(k+nseln)+1] = -Hm[k]*t.y;
						N[3*(k+nseln)+2] = -Hm[k]*t.z;
					}
					
					for (k=0; k<ndof; ++k) fe[k] += N[k]*detJ[j]*w[j];
					
					// assemble the global residual
					R.Assemble(en, LM, fe);
					
					// do the biphasic stuff
					// TODO: I should only do this when the node is actually in contact
					if (sporo && mporo && pt.m_pme)
					{
						// calculate nr of pressure dofs
						int ndof = nseln + nmeln;
						
						// calculate the flow rate
						double epsp = m_epsp*pt.m_epsp;
						
						double wn = pt.m_Lmp + epsp*pt.m_pg;
						
						// fill the LM
						LM.resize(ndof);
						for (k=0; k<nseln; ++k) LM[k        ] = sLM[3*nseln+k];
						for (k=0; k<nmeln; ++k) LM[k + nseln] = mLM[3*nmeln+k];
						
						// fill the force array
						fe.resize(ndof);
						zero(fe);
						for (k=0; k<nseln; ++k) N[k      ] =  Hs[k];
						for (k=0; k<nmeln; ++k) N[k+nseln] = -Hm[k];
						
						for (k=0; k<ndof; ++k) fe[k] += dt*wn*N[k]*detJ[j]*w[j];
						
						// assemble residual
						R.Assemble(en, LM, fe);
					}
				}
			}
		}
	}
}

//-----------------------------------------------------------------------------
void FETiedBiphasicInterface::StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp)
{
	int i, j, k, l;
	vector<int> sLM, mLM, LM, en;
	const int MN = FEElement::MAX_NODES;
	double detJ[MN], w[MN], *Hs, Hm[MN], pt[MN], dpr[MN], dps[MN];
	FEElementMatrix ke;

	// get time step
    double dt = GetFEModel()->GetTime().timeIncrement;
	
	// get the mesh
	FEMesh* pm = m_ss.GetMesh();
	
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
		FETiedBiphasicSurface& ss = (np == 0? m_ss : m_ms);
		FETiedBiphasicSurface& ms = (np == 0? m_ms : m_ss);
		
		// loop over all primary surface elements
		for (i=0; i<ss.Elements(); ++i)
		{
			// get the next element
			FESurfaceElement& se = ss.Element(i);
			
			bool sporo = ss.m_poro[i];
			
			// get nr of nodes and integration points
			int nseln = se.Nodes();
			int nint = se.GaussPoints();
			
			// nodal pressures
			double pn[FEElement::MAX_NODES];
			for (j=0; j<nseln; ++j) pn[j] = ss.GetMesh()->Node(se.m_node[j]).get(m_dofP);
			
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
			}
			
			// loop over all integration points
			for (j=0; j<nint; ++j)
			{
				FETiedBiphasicSurface::Data& pt = static_cast<FETiedBiphasicSurface::Data&>(*se.GetMaterialPoint(j));

				// get the secondary surface element
				FESurfaceElement* pme = pt.m_pme;
				if (pme)
				{
					FESurfaceElement& me = *pme;
					
					bool mporo = ms.m_poro[pme->m_lid];
					
					// get the nr of secondary surface nodes
					int nmeln = me.Nodes();
					
					// nodal pressure
					double pm[FEElement::MAX_NODES];
					if (mporo) for (k=0; k<nmeln; ++k) pm[k] = ms.GetMesh()->Node(me.m_node[k]).get(m_dofP);
					
					// copy the LM vector
					ms.UnpackLM(me, mLM);
					
					int ndpn;	// number of dofs per node
					int ndof;	// number of dofs in stiffness matrix
					
					if (sporo && mporo) {
						// calculate degrees of freedom for biphasic-on-biphasic contact
						ndpn = 4;
						ndof = ndpn*(nseln+nmeln);
						
						// build the LM vector
						LM.resize(ndof);
						
						for (k=0; k<nseln; ++k)
						{
							LM[4*k  ] = sLM[3*k  ];			// x-dof
							LM[4*k+1] = sLM[3*k+1];			// y-dof
							LM[4*k+2] = sLM[3*k+2];			// z-dof
							LM[4*k+3] = sLM[3*nseln+k];		// p-dof
						}
						for (k=0; k<nmeln; ++k)
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
					vec3d dg = pt.m_dg;
					
					// lagrange multiplier
					vec3d Lm = pt.m_Lmd;
					
					// penalty 
					double eps = m_epsn*pt.m_epsn;
					
					// contact traction
					vec3d t = Lm + dg*eps;
					
					// create the stiffness matrix
					ke.resize(ndof, ndof); ke.zero();
					
					// --- S O L I D - S O L I D   C O N T A C T ---
					
					// a. I-term
					//------------------------------------
					
					for (k=0; k<nseln; ++k) {
						for (l=0; l<nseln; ++l)
						{
							ke[ndpn*k    ][ndpn*l    ] += eps*Hs[k]*Hs[l]*detJ[j]*w[j];
							ke[ndpn*k + 1][ndpn*l + 1] += eps*Hs[k]*Hs[l]*detJ[j]*w[j];
							ke[ndpn*k + 2][ndpn*l + 2] += eps*Hs[k]*Hs[l]*detJ[j]*w[j];
						}
						for (l=0; l<nmeln; ++l)
						{
							ke[ndpn*k    ][ndpn*(nseln+l)    ] += -eps*Hs[k]*Hm[l]*detJ[j]*w[j];
							ke[ndpn*k + 1][ndpn*(nseln+l) + 1] += -eps*Hs[k]*Hm[l]*detJ[j]*w[j];
							ke[ndpn*k + 2][ndpn*(nseln+l) + 2] += -eps*Hs[k]*Hm[l]*detJ[j]*w[j];
						}
					}
					
					for (k=0; k<nmeln; ++k) {
						for (l=0; l<nseln; ++l)
						{
							ke[ndpn*(nseln+k)    ][ndpn*l    ] += -eps*Hm[k]*Hs[l]*detJ[j]*w[j];
							ke[ndpn*(nseln+k) + 1][ndpn*l + 1] += -eps*Hm[k]*Hs[l]*detJ[j]*w[j];
							ke[ndpn*(nseln+k) + 2][ndpn*l + 2] += -eps*Hm[k]*Hs[l]*detJ[j]*w[j];
						}
						for (l=0; l<nmeln; ++l)
						{
							ke[ndpn*(nseln+k)    ][ndpn*(nseln+l)    ] += eps*Hm[k]*Hm[l]*detJ[j]*w[j];
							ke[ndpn*(nseln+k) + 1][ndpn*(nseln+l) + 1] += eps*Hm[k]*Hm[l]*detJ[j]*w[j];
							ke[ndpn*(nseln+k) + 2][ndpn*(nseln+l) + 2] += eps*Hm[k]*Hm[l]*detJ[j]*w[j];
						}
					}
					
					// b. A-term
					//-------------------------------------
					
					double* Gr = se.Gr(j);
					double* Gs = se.Gs(j);
					vec3d gs[2];
					ss.CoBaseVectors(se, j, gs);
					
					vec3d as[FEElement::MAX_NODES];
					mat3d As[FEElement::MAX_NODES];
					for (l=0; l<nseln; ++l) {
						as[l] = nu ^ (gs[1]*Gr[l] - gs[0]*Gs[l]);
						As[l] = t & as[l];
					}
					
					if (!m_bsymm)
					{
						// non-symmetric
						for (k=0; k<nseln; ++k) {
							for (l=0; l<nseln; ++l)
							{
								ke[ndpn*k    ][ndpn*l    ] += Hs[k]*As[l](0,0)*w[j];
								ke[ndpn*k    ][ndpn*l + 1] += Hs[k]*As[l](0,1)*w[j];
								ke[ndpn*k    ][ndpn*l + 2] += Hs[k]*As[l](0,2)*w[j];

								ke[ndpn*k + 1][ndpn*l    ] += Hs[k]*As[l](1,0)*w[j];
								ke[ndpn*k + 1][ndpn*l + 1] += Hs[k]*As[l](1,1)*w[j];
								ke[ndpn*k + 1][ndpn*l + 2] += Hs[k]*As[l](1,2)*w[j];

								ke[ndpn*k + 2][ndpn*l    ] += Hs[k]*As[l](2,0)*w[j];
								ke[ndpn*k + 2][ndpn*l + 1] += Hs[k]*As[l](2,1)*w[j];
								ke[ndpn*k + 2][ndpn*l + 2] += Hs[k]*As[l](2,2)*w[j];
							}
						}
						
						for (k=0; k<nmeln; ++k) {
							for (l=0; l<nseln; ++l)
							{
								ke[ndpn*(nseln+k)    ][ndpn*l    ] += -Hm[k]*As[l](0,0)*w[j];
								ke[ndpn*(nseln+k)    ][ndpn*l + 1] += -Hm[k]*As[l](0,1)*w[j];
								ke[ndpn*(nseln+k)    ][ndpn*l + 2] += -Hm[k]*As[l](0,2)*w[j];

								ke[ndpn*(nseln+k) + 1][ndpn*l    ] += -Hm[k]*As[l](1,0)*w[j];
								ke[ndpn*(nseln+k) + 1][ndpn*l + 1] += -Hm[k]*As[l](1,1)*w[j];
								ke[ndpn*(nseln+k) + 1][ndpn*l + 2] += -Hm[k]*As[l](1,2)*w[j];

								ke[ndpn*(nseln+k) + 2][ndpn*l    ] += -Hm[k]*As[l](2,0)*w[j];
								ke[ndpn*(nseln+k) + 2][ndpn*l + 1] += -Hm[k]*As[l](2,1)*w[j];
								ke[ndpn*(nseln+k) + 2][ndpn*l + 2] += -Hm[k]*As[l](2,2)*w[j];
							}
						}
						
					}
					else 
					{
						// symmetric
						for (k=0; k<nseln; ++k) {
							for (l=0; l<nseln; ++l)
							{
								ke[ndpn*k    ][ndpn*l    ] += 0.5*(Hs[k]*As[l](0,0)+Hs[l]*As[k](0,0))*w[j];
								ke[ndpn*k    ][ndpn*l + 1] += 0.5*(Hs[k]*As[l](0,1)+Hs[l]*As[k](1,0))*w[j];
								ke[ndpn*k    ][ndpn*l + 2] += 0.5*(Hs[k]*As[l](0,2)+Hs[l]*As[k](2,0))*w[j];
								
								ke[ndpn*k + 1][ndpn*l    ] += 0.5*(Hs[k]*As[l](1,0)+Hs[l]*As[k](0,1))*w[j];
								ke[ndpn*k + 1][ndpn*l + 1] += 0.5*(Hs[k]*As[l](1,1)+Hs[l]*As[k](1,1))*w[j];
								ke[ndpn*k + 1][ndpn*l + 2] += 0.5*(Hs[k]*As[l](1,2)+Hs[l]*As[k](2,1))*w[j];
								
								ke[ndpn*k + 2][ndpn*l    ] += 0.5*(Hs[k]*As[l](2,0)+Hs[l]*As[k](0,2))*w[j];
								ke[ndpn*k + 2][ndpn*l + 1] += 0.5*(Hs[k]*As[l](2,1)+Hs[l]*As[k](1,2))*w[j];
								ke[ndpn*k + 2][ndpn*l + 2] += 0.5*(Hs[k]*As[l](2,2)+Hs[l]*As[k](2,2))*w[j];
							}
						}
						
						for (k=0; k<nmeln; ++k) {
							for (l=0; l<nseln; ++l)
							{
								ke[ndpn*(nseln+k)    ][ndpn*l    ] += -0.5*Hm[k]*As[l](0,0)*w[j];
								ke[ndpn*(nseln+k)    ][ndpn*l + 1] += -0.5*Hm[k]*As[l](0,1)*w[j];
								ke[ndpn*(nseln+k)    ][ndpn*l + 2] += -0.5*Hm[k]*As[l](0,2)*w[j];
								
								ke[ndpn*(nseln+k) + 1][ndpn*l    ] += -0.5*Hm[k]*As[l](1,0)*w[j];
								ke[ndpn*(nseln+k) + 1][ndpn*l + 1] += -0.5*Hm[k]*As[l](1,1)*w[j];
								ke[ndpn*(nseln+k) + 1][ndpn*l + 2] += -0.5*Hm[k]*As[l](1,2)*w[j];
								
								ke[ndpn*(nseln+k) + 2][ndpn*l    ] += -0.5*Hm[k]*As[l](2,0)*w[j];
								ke[ndpn*(nseln+k) + 2][ndpn*l + 1] += -0.5*Hm[k]*As[l](2,1)*w[j];
								ke[ndpn*(nseln+k) + 2][ndpn*l + 2] += -0.5*Hm[k]*As[l](2,2)*w[j];
							}
						}
						
						for (k=0; k<nseln; ++k) {
							for (l=0; l<nmeln; ++l)
							{
								ke[ndpn*k    ][ndpn*(nseln+l)    ] += -0.5*Hm[l]*As[k](0,0)*w[j];
								ke[ndpn*k    ][ndpn*(nseln+l) + 1] += -0.5*Hm[l]*As[k](1,0)*w[j];
								ke[ndpn*k    ][ndpn*(nseln+l) + 2] += -0.5*Hm[l]*As[k](2,0)*w[j];
								
								ke[ndpn*k + 1][ndpn*(nseln+l)    ] += -0.5*Hm[l]*As[k](0,1)*w[j];
								ke[ndpn*k + 1][ndpn*(nseln+l) + 1] += -0.5*Hm[l]*As[k](1,1)*w[j];
								ke[ndpn*k + 1][ndpn*(nseln+l) + 2] += -0.5*Hm[l]*As[k](2,1)*w[j];
								
								ke[ndpn*k + 2][ndpn*(nseln+l)    ] += -0.5*Hm[l]*As[k](0,2)*w[j];
								ke[ndpn*k + 2][ndpn*(nseln+l) + 1] += -0.5*Hm[l]*As[k](1,2)*w[j];
								ke[ndpn*k + 2][ndpn*(nseln+l) + 2] += -0.5*Hm[l]*As[k](2,2)*w[j];
							}
						}
					}

					
					// --- B I P H A S I C   S T I F F N E S S ---
					if (sporo && mporo)
					{
						double epsp = (pt.m_pme) ? m_epsp*pt.m_epsp : 0.;
						
						// --- S O L I D - P R E S S U R E   C O N T A C T ---
						
						// b. A-term
						//-------------------------------------

						double wn = pt.m_Lmp + epsp*pt.m_pg;
						
						if (!m_bsymm)
						{
							// non-symmetric
							for (k=0; k<nseln; ++k)
								for (l=0; l<nseln; ++l) {
								{
									ke[4*k + 3][4*l  ] += dt*w[j]*wn*Hs[k]*as[l].x;
									ke[4*k + 3][4*l+1] += dt*w[j]*wn*Hs[k]*as[l].y;
									ke[4*k + 3][4*l+2] += dt*w[j]*wn*Hs[k]*as[l].z;
								}
							}
							for (k=0; k<nmeln; ++k)
								for (l=0; l<nseln; ++l) {
									{
										ke[4*(k+nseln) + 3][4*l  ] += -dt*w[j]*wn*Hm[k]*as[l].x;
										ke[4*(k+nseln) + 3][4*l+1] += -dt*w[j]*wn*Hm[k]*as[l].y;
										ke[4*(k+nseln) + 3][4*l+2] += -dt*w[j]*wn*Hm[k]*as[l].z;
									}
								}
						}
						else 
						{
							// symmetric
							for (k=0; k<nseln; ++k)
								for (l=0; l<nseln; ++l) {
									{
										ke[4*k + 3][4*l  ] += dt*w[j]*wn*0.5*(Hs[k]*as[l].x+Hs[l]*as[k].x);
										ke[4*k + 3][4*l+1] += dt*w[j]*wn*0.5*(Hs[k]*as[l].y+Hs[l]*as[k].y);
										ke[4*k + 3][4*l+2] += dt*w[j]*wn*0.5*(Hs[k]*as[l].z+Hs[l]*as[k].z);
									}
								}
							for (k=0; k<nmeln; ++k)
								for (l=0; l<nseln; ++l) {
									{
										ke[4*(k+nseln) + 3][4*l  ] += -dt*w[j]*wn*0.5*Hm[k]*as[l].x;
										ke[4*(k+nseln) + 3][4*l+1] += -dt*w[j]*wn*0.5*Hm[k]*as[l].y;
										ke[4*(k+nseln) + 3][4*l+2] += -dt*w[j]*wn*0.5*Hm[k]*as[l].z;
									}
								}
							for (k=0; k<nseln; ++k)
								for (l=0; l<nmeln; ++l) {
									{
										ke[4*k + 3][4*(nseln+l)  ] += -dt*w[j]*wn*0.5*Hm[l]*as[k].x;
										ke[4*k + 3][4*(nseln+l)+1] += -dt*w[j]*wn*0.5*Hm[l]*as[k].y;
										ke[4*k + 3][4*(nseln+l)+2] += -dt*w[j]*wn*0.5*Hm[l]*as[k].z;
									}
								}
						}

						
						// --- P R E S S U R E - P R E S S U R E   C O N T A C T ---
						
						for (k=0; k<nseln; ++k) {
							for (l=0; l<nseln; ++l)
								ke[4*k + 3][4*l+3] += -dt*epsp*w[j]*detJ[j]*Hs[k]*Hs[l];
							for (l=0; l<nmeln; ++l)
								ke[4*k + 3][4*(nseln+l)+3] += dt*epsp*w[j]*detJ[j]*Hs[k]*Hm[l];
						}
						
						for (k=0; k<nmeln; ++k) {
							for (l=0; l<nseln; ++l)
								ke[4*(nseln+k)+3][4*l + 3] += dt*epsp*w[j]*detJ[j]*Hm[k]*Hs[l];
							for (l=0; l<nmeln; ++l)
								ke[4*(nseln+k)+3][4*(nseln+l) + 3] += -dt*epsp*w[j]*detJ[j]*Hm[k]*Hm[l];
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
bool FETiedBiphasicInterface::Augment(int naug, const FETimeInfo& tp)
{
	// make sure we need to augment
	if (m_laugon != 1) return true;

	int i;
	vec3d Ln;
	double Lp;
	bool bconv = true;
	
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
			FETiedBiphasicSurface::Data& ds = static_cast<FETiedBiphasicSurface::Data&>(*el.GetMaterialPoint(j));
			normL0 += ds.m_Lmd*ds.m_Lmd;
        }
    }
    for (int i=0; i<NM; ++i)
    {
		FESurfaceElement& el = m_ms.Element(i);
		for (int j = 0; j<el.GaussPoints(); ++j)
		{
			FETiedBiphasicSurface::Data& dm = static_cast<FETiedBiphasicSurface::Data&>(*el.GetMaterialPoint(j));
			normL0 += dm.m_Lmd*dm.m_Lmd;
        }
    }
	
	// b. gap component
	// (is calculated during update)
	double maxgap = 0;
	double maxpg = 0;
	
	// update Lagrange multipliers
	double normL1 = 0, eps, epsp;
	for (int i = 0; i<NS; ++i)
	{
		FESurfaceElement& el = m_ss.Element(i);
		for (int j = 0; j<el.GaussPoints(); ++j)
		{
			FETiedBiphasicSurface::Data& ds = static_cast<FETiedBiphasicSurface::Data&>(*el.GetMaterialPoint(j));

            // update Lagrange multipliers on primary surface
            eps = m_epsn*ds.m_epsn;
            ds.m_Lmd = ds.m_Lmd + ds.m_dg*eps;
            
            normL1 += ds.m_Lmd*ds.m_Lmd;
            
            if (m_ss.m_bporo) {
                Lp = 0;
                if (ds.m_pme) {
                    epsp = m_epsp*ds.m_epsp;
                    Lp = ds.m_Lmp + epsp*ds.m_pg;
                    maxpg = max(maxpg,fabs(ds.m_pg));
                    normDP += ds.m_pg*ds.m_pg;
                }
                ds.m_Lmp = Lp;
            }
            
            maxgap = max(maxgap,sqrt(ds.m_dg*ds.m_dg));
        }
    }
	
	for (i=0; i<NM; ++i)
	{
		FESurfaceElement& el = m_ms.Element(i);
		for (int j = 0; j<el.GaussPoints(); ++j)
		{
			FETiedBiphasicSurface::Data& dm = static_cast<FETiedBiphasicSurface::Data&>(*el.GetMaterialPoint(j));

            // update Lagrange multipliers on secondary surface
            eps = m_epsn*dm.m_epsn;
            dm.m_Lmd = dm.m_Lmd + dm.m_dg*eps;
            
            normL1 += dm.m_Lmd*dm.m_Lmd;
            
            if (m_ms.m_bporo) {
                Lp = 0;
                if (dm.m_pme) {
                    epsp = m_epsp*dm.m_epsp;
                    Lp = dm.m_Lmp + epsp*dm.m_pg;
                    maxpg = max(maxpg,fabs(dm.m_pg));
                    normDP += dm.m_pg*dm.m_pg;
                }
                dm.m_Lmp = Lp;
            }
            
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
	if ((m_ptol > 0) && (bporo && maxpg > m_ptol)) bconv = false;
	
	if ((m_atol > 0) && (lnorm > m_atol)) bconv = false;
	if ((m_atol > 0) && (pnorm > m_atol)) bconv = false;

	if (naug < m_naugmin ) bconv = false;
	if (naug >= m_naugmax) bconv = true;

	feLog(" tied biphasic interface # %d\n", GetID());
	feLog("                        CURRENT        REQUIRED\n");
	feLog("    D multiplier : %15le", lnorm); if (m_atol > 0) feLog("%15le\n", m_atol); else feLog("       ***\n");
	if (bporo) { feLog("    P gap        : %15le", pnorm); if (m_atol > 0) feLog("%15le\n", m_atol); else feLog("       ***\n"); }
	
	feLog("    maximum gap  : %15le", maxgap);
	if (m_gtol > 0) feLog("%15le\n", m_gtol); else feLog("       ***\n");
	if (bporo) {
		feLog("    maximum pgap : %15le", maxpg);
		if (m_ptol > 0) feLog("%15le\n", m_ptol); else feLog("       ***\n");
	}
	
	return bconv;
}

//-----------------------------------------------------------------------------
void FETiedBiphasicInterface::Serialize(DumpStream &ar)
{
	// store contact data
	FEContactInterface::Serialize(ar);

	// store contact surface data
	m_ms.Serialize(ar);
	m_ss.Serialize(ar);

	// serialize pointers
	SerializeElementPointers(m_ss, m_ms, ar);
	SerializeElementPointers(m_ms, m_ss, ar);
}
