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
#include "FESlidingInterface.h"
#include <FECore/FEShellDomain.h>
#include "FECore/FEClosestPointProjection.h"
#include "FECore/FEModel.h"
#include "FECore/FEGlobalMatrix.h"
#include "FECore/log.h"
#include <FECore/FELinearSystem.h>
#include <FECore/FEAnalysis.h>

FESlidingSurface::FESlidingPoint::FESlidingPoint()
{
	m_gap = 0.0;
	m_nu = vec3d(0,0,0);
	m_rs = vec2d(0,0);
	m_rsp = vec2d(0,0);
	m_Lm = 0.0;
	m_M.zero();
	m_Lt = vec2d(0,0);
	m_off = 0.0;
	m_eps = 1.0;
	m_Ln = 0.0;
}

void FESlidingSurface::FESlidingPoint::Serialize(DumpStream& ar)
{
	FEContactMaterialPoint::Serialize(ar);

	ar & m_nu & m_eps & m_off;
	ar & m_rs & m_rsp;
	ar & m_Lm & m_Lt;
	ar & m_M;
}

void FESlidingSurface::FESlidingPoint::Init()
{
	FEContactMaterialPoint::Init();
	m_gap = 0.0;
	m_nu = vec3d(0, 0, 0);
	m_rs = vec2d(0, 0);
	m_rsp = vec2d(0, 0);
	m_Lm = 0.0;
	m_M.zero();
	m_Lt = vec2d(0, 0);
	m_off = 0.0;
	m_eps = 1.0;
}

//-----------------------------------------------------------------------------
// Define sliding interface parameters
BEGIN_FECORE_CLASS(FESlidingInterface, FEContactInterface)
	ADD_PARAMETER(m_atol         , "tolerance"    );
	ADD_PARAMETER(m_eps          , "penalty"      );
	ADD_PARAMETER(m_bautopen     , "auto_penalty" );
	ADD_PARAMETER(m_btwo_pass    , "two_pass"     );
	ADD_PARAMETER(m_gtol         , "gaptol"       );
	ADD_PARAMETER(m_mu           , "fric_coeff"   );
	ADD_PARAMETER(m_epsf         , "fric_penalty" );
	ADD_PARAMETER(m_naugmin      , "minaug"       );
	ADD_PARAMETER(m_naugmax      , "maxaug"       );
	ADD_PARAMETER(m_stol         , "search_tol"   );
	ADD_PARAMETER(m_ktmult       , "ktmult"       );
	ADD_PARAMETER(m_knmult       , "knmult"       );
	ADD_PARAMETER(m_breloc       , "node_reloc"   );
	ADD_PARAMETER(m_nsegup       , "seg_up"       );
	ADD_PARAMETER(m_sradius      , "search_radius");
	ADD_PARAMETER(m_bupdtpen     , "update_penalty");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FESlidingSurface::FESlidingSurface(FEModel* pfem) : FEContactSurface(pfem) {}

//-----------------------------------------------------------------------------
//! build the matrix profile for use in the stiffness matrix
void FESlidingInterface::BuildMatrixProfile(FEGlobalMatrix& K)
{
	// TODO: this is currently for max 6 nodes (hence 7=6+1)
	vector<int> lm(6*7);

	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// get the DOFS
	const int dof_X = fem.GetDOFIndex("x");
	const int dof_Y = fem.GetDOFIndex("y");
	const int dof_Z = fem.GetDOFIndex("z");
	const int dof_RU = fem.GetDOFIndex("Ru");
	const int dof_RV = fem.GetDOFIndex("Rv");
	const int dof_RW = fem.GetDOFIndex("Rw");

	int npass = (m_btwo_pass?2:1);
	for (int np=0; np<npass; ++np)
	{
		FESlidingSurface& ss = (np==0? m_ss : m_ms);
		FESlidingSurface& ms = (np==0? m_ms : m_ss);

		for (int j=0; j<ss.Nodes(); ++j)
		{
			FESurfaceElement* pe = ss.m_data[j].m_pme;

			if (pe != 0)
			{
				FESurfaceElement& me = *pe;
				int* en = &me.m_node[0];

				// Note that we need to grab the rigid degrees of freedom as well
				// this is in case one of the nodes belongs to a rigid body.
				int n = me.Nodes();
				if (n == 3)
				{
					lm[6*(3+1)  ] = -1;lm[6*(3+2)  ] = -1;lm[6*(3+3)  ] = -1;
					lm[6*(3+1)+1] = -1;lm[6*(3+2)+1] = -1;lm[6*(3+3)+1] = -1;
					lm[6*(3+1)+2] = -1;lm[6*(3+2)+2] = -1;lm[6*(3+3)+2] = -1;
					lm[6*(3+1)+3] = -1;lm[6*(3+2)+3] = -1;lm[6*(3+3)+3] = -1;
					lm[6*(3+1)+4] = -1;lm[6*(3+2)+4] = -1;lm[6*(3+3)+4] = -1;
					lm[6*(3+1)+5] = -1;lm[6*(3+2)+5] = -1;lm[6*(3+3)+5] = -1;
				}
				if (n == 4)
				{
					lm[6*(4+1)  ] = -1;lm[6*(4+2)  ] = -1;
					lm[6*(4+1)+1] = -1;lm[6*(4+2)+1] = -1;
					lm[6*(4+1)+2] = -1;lm[6*(4+2)+2] = -1;
					lm[6*(4+1)+3] = -1;lm[6*(4+2)+3] = -1;
					lm[6*(4+1)+4] = -1;lm[6*(4+2)+4] = -1;
					lm[6*(4+1)+5] = -1;lm[6*(4+2)+5] = -1;
				}

				lm[0] = ss.Node(j).m_ID[dof_X];
				lm[1] = ss.Node(j).m_ID[dof_Y];
				lm[2] = ss.Node(j).m_ID[dof_Z];
				lm[3] = ss.Node(j).m_ID[dof_RU];
				lm[4] = ss.Node(j).m_ID[dof_RV];
				lm[5] = ss.Node(j).m_ID[dof_RW];

				for (int k=0; k<n; ++k)
				{
					vector<int>& id = mesh.Node(en[k]).m_ID;
					lm[6*(k+1)  ] = id[dof_X];
					lm[6*(k+1)+1] = id[dof_Y];
					lm[6*(k+1)+2] = id[dof_Z];
					lm[6*(k+1)+3] = id[dof_RU];
					lm[6*(k+1)+4] = id[dof_RV];
					lm[6*(k+1)+5] = id[dof_RW];
				}

				K.build_add(lm);
			}
		}
	}
}

//-----------------------------------------------------------------------------
//! Creates a surface for use with a sliding interface. All surface data
//! structures are allocated.
//! Note that it is assumed that the element array is already created
//! and initialized.

bool FESlidingSurface::Init()
{
	// always intialize base class first!
	if (FEContactSurface::Init() == false) return false;

	// make sure the sibling surface has been set
	assert(m_pSibling);

	// get the number of nodes
	int nn = Nodes();

	// allocate integration point data
	m_data.resize(nn);
	for (int i = 0; i < nn; ++i)
	{
		FESlidingPoint& d = m_data[i];
		d.Init();
	}

	// we calculate the gap offset values
	// This value is used to take the shell thickness into account
	// note that we force rigid shells to have zero thickness
	FEMesh& m = *m_pMesh;
	vector<double> tag(m.Nodes());
	zero(tag);
	for (int nd=0; nd<m.Domains(); ++nd)
	{
		FEShellDomain* psd = dynamic_cast<FEShellDomain*>(&m.Domain(nd));
		if (psd)
		{
			for (int i=0; i<psd->Elements(); ++i)
			{
				FEShellElement& el = psd->Element(i);
				int n = el.Nodes();
				for (int j=0; j<n; ++j) tag[el.m_node[j]] = 0.5*el.m_h0[j];
			}
		}
	}
	for (int i=0; i<nn; ++i) m_data[i].m_off = tag[NodeIndex(i)];

	return true;
}

//-----------------------------------------------------------------------------
//! 
vec3d FESlidingSurface::traction(int inode)
{
	vec3d t(0,0,0);
	if (m_data[inode].m_pme)
	{
		FESurfaceElement& el = *m_data[inode].m_pme;
		double Tn = m_data[inode].m_Lm;
		double T1 = -m_data[inode].m_Lt[0];
		double T2 = -m_data[inode].m_Lt[1];
		double r = m_data[inode].m_rs[0];
		double s = m_data[inode].m_rs[1];
        
		vec3d tn = m_data[inode].m_nu*Tn, tt;
		vec3d e[2];
		ContraBaseVectors(el, r, s, e);
		tt = e[0]*T1 + e[1]*T2;
		t = tn + tt;
	}
    
	return t;
}

//-----------------------------------------------------------------------------
vec3d FESlidingSurface::GetContactForce()
{
	const int MN = FEElement::MAX_NODES;
	double Tn[MN],T1[MN],T2[MN];
	
	// initialize contact force
	vec3d f(0,0,0);
	
	// loop over all elements of the surface
	for (int n=0; n<Elements(); ++n)
	{
		FESurfaceElement& el = Element(n);
		int nseln = el.Nodes();
		
		// nodal contact pressures and frictional tractions
		for (int i=0; i<nseln; ++i) {
            Tn[i] =  m_data[el.m_lnode[i]].m_Ln;
            T1[i] = -m_data[el.m_lnode[i]].m_Lt[0];
            T2[i] = -m_data[el.m_lnode[i]].m_Lt[1];
        }
		int nint = el.GaussPoints();
		
		// evaluate the contact force for that element
		for (int i=0; i<nint; ++i)
		{
			// area in reference configuration
			vec3d g0[2],g[2];
			double r = el.gr(i);
			double s = el.gs(i);
			CoBaseVectors0(el, r, s, g0);
			double A = (g0[0] ^ g0[1]).unit();
			// traction components at integration point
            double t1 = el.eval(T1,i);
            double t2 = el.eval(T2,i);
			double t3 = el.eval(Tn,i);
			// unit normal vector
			vec3d n = SurfaceNormal(el, i);
            // contravariant basis in spatial frame
            ContraBaseVectors(el, r, s, g);
            // Piola traction
            vec3d t = g[0]*t1 + g[1]*t2 + n*t3;
			// gauss weight
			double w = el.GaussWeights()[i];
			// contact force
			f += t*(w*A);
		}
	}
	
	return f;
}

//-----------------------------------------------------------------------------
double FESlidingSurface::GetContactArea()
{
	const int MN = FEElement::MAX_NODES;
	double Tn[MN];
    
	// initialize contact area
	double a = 0;
	
	// loop over all elements of the primary surface
	for (int n=0; n<Elements(); ++n)
	{
		FESurfaceElement& el = Element(n);
		int nint = el.GaussPoints();
		
		int nseln = el.Nodes();
		
		// nodal contact pressures
		for (int i=0; i<nseln; ++i) {
            Tn[i] = m_data[el.m_lnode[i]].m_Ln;
        }
        
		// evaluate the contact force for that element
		for (int i=0; i<nint; ++i)
		{
			// get data for this integration point
			double Ln = el.eval(Tn,i);
            double s = (Ln > 0) ? 1 : 0;
            
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
void FESlidingSurface::Serialize(DumpStream& ar)
{
	// serialize base class
	FEContactSurface::Serialize(ar);

	// serialize data
	ar & m_data;
}

//-----------------------------------------------------------------------------
void FESlidingSurface::GetContactTraction(int nface, vec3d& pt)
{
    FESurfaceElement& el = Element(nface);
    int ne = el.Nodes();
    pt = vec3d(0,0,0);
    for (int k=0; k<ne; ++k)
    {
        int nj = el.m_lnode[k];
        if (m_data[nj].m_gap > 0) pt += m_data[nj].m_nu*m_data[nj].m_Ln;
    }
    pt /= ne;
}

//-----------------------------------------------------------------------------
void FESlidingSurface::GetNodalContactPressure(int nface, double* pg)
{
	FESurfaceElement& f = Element(nface);
	int ne = f.Nodes();
	for (int j=0; j<ne; ++j) pg[j] = m_data[f.m_lnode[j]].m_Ln;
}

//-----------------------------------------------------------------------------
void FESlidingSurface::GetNodalContactTraction(int nface, vec3d* tn)
{
	FESurfaceElement& e = Element(nface);
	int ne = e.Nodes();
	for (int j=0; j<ne; ++j)
	{
		int nj = e.m_lnode[j];
		double gi = m_data[nj].m_gap;
		double Li = m_data[nj].m_Ln;
		vec3d ti  = m_data[nj].m_nu;
		if (gi > 0) tn[j] = ti*Li; else tn[j] = vec3d(0,0,0);
	}
}

///////////////////////////////////////////////////////////////////////////////
// FESlidingInterface
///////////////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
//! constructor
FESlidingInterface::FESlidingInterface(FEModel* pfem) : FEContactInterface(pfem), m_ss(pfem), m_ms(pfem)
{
	static int count = 1;
	SetID(count++);

	m_mu = 0;
	m_epsf = 0;

	m_naugmin = 0;
	m_naugmax = 10;

	m_gtol = 0;

	m_stol = 0.01;

	m_ktmult = 0;
	m_knmult = 1;

	m_breloc = false;

	m_bupdtpen = false;

	m_nsegup = 0;	// always do segment updates
	m_bautopen = false;	// don't use auto-penalty
	m_btwo_pass = false; // don't use two-pass
	m_sradius = 0;				// no search radius limitation

	// set parents
	m_ms.SetContactInterface(this);
	m_ss.SetContactInterface(this);

	// set the siblings
	m_ms.SetSibling(&m_ss);
	m_ss.SetSibling(&m_ms);
};

//-----------------------------------------------------------------------------
//! Calculates the auto penalty factor

void FESlidingInterface::CalcAutoPenalty(FESlidingSurface& s)
{
	// zero penalty values
	for (int i=0; i<s.Nodes(); ++i) s.m_data[i].m_eps = 0.0;

	// get the mesh
	FEMesh& mesh = *s.GetMesh();

	// get the node element list for this surface
	FENodeElemList NEL;
	NEL.Create(s);

	// loop over all surface elements
	for (int i=0; i<s.Elements(); ++i)
	{
		// get the next face
		FESurfaceElement& face = s.Element(i);

		// we need a measure for the modulus
		double eps = AutoPenalty(face, s);

		// distribute values over nodes
		for (int k=0; k<face.Nodes(); ++k)
		{
			int m = face.m_lnode[k];
			s.m_data[m].m_eps += eps;
		}
	}

	// scale values according to valence (TODO: Why are we doing this?)
	for (int i=0; i<s.Nodes(); ++i) s.m_data[i].m_eps /= NEL.Valence(i);
}

//-----------------------------------------------------------------------------
//! Initializes the sliding interface data

bool FESlidingInterface::Init()
{
	// set data
	m_bfirst = true;
	m_normg0 = 0.0;

	// create the surfaces
	if (m_ss.Init() == false) return false;
	if (m_ms.Init() == false) return false;

	return true;
}

//-----------------------------------------------------------------------------
void FESlidingInterface::Activate()
{
	// don't forget to call the base class
	FEContactInterface::Activate();

	// project primary surface onto secondary surface
	ProjectSurface(m_ss, m_ms, true, m_breloc);
	if (m_bautopen) CalcAutoPenalty(m_ss);

	// for two-pass algorithms we repeat the previous
	// two steps with primary and secondary surface switched
	if (m_btwo_pass)
	{
		ProjectSurface(m_ms, m_ss, true);
		if (m_bautopen) CalcAutoPenalty(m_ms);
	}
}

//-----------------------------------------------------------------------------
//!  Projects the primary surface onto the secondary surface.
//!  That is, for each primary surface node we determine the closest
//!  secondary surface element and the projection of that node onto
//!  this element.

//! \todo this function needs to identify the different types of contact:
//!   1/ first contact
//!   2/ crossing of element boundary
//!	  3/ contact termination 
//!			either by failure to find projection or when g < tolerance

void FESlidingInterface::ProjectSurface(FESlidingSurface& ss, FESlidingSurface& ms, bool bupseg, bool bmove)
{
	// node projection data
	double r, s;
	vec3d q;

	FEClosestPointProjection cpp(ms);
	cpp.SetTolerance(m_stol);
	cpp.SetSearchRadius(m_sradius);
	cpp.HandleSpecialCases(true);
	cpp.Init();

	// loop over all primary surface nodes
	for (int i=0; i<ss.Nodes(); ++i)
	{
		// get the node
		FENode& node = ss.Node(i);

		// get the nodal position
		vec3d x = node.m_rt;

		// get the global node number
		int m = ss.NodeIndex(i);

		// get the previous secondary surface element (if any)
		FESurfaceElement* pme = ss.m_data[i].m_pme;

		// If the node is in contact, let's see if the node still is 
		// on the same element
		if (pme != 0)
		{
			FESurfaceElement& mel = *pme;

			r = ss.m_data[i].m_rs[0];
			s = ss.m_data[i].m_rs[1];

			q = ms.ProjectToSurface(mel, x, r, s);
			ss.m_data[i].m_rs[0] = r;
			ss.m_data[i].m_rs[1] = s;

			// we only check when we can update the segments
			// otherwise, we just stick with this element, even
			// if the node is no longer inside it.
			if (bupseg)
			{
				if (!ms.IsInsideElement(mel, r, s, m_stol))
				{
					// see if the node might have moved to another element
					FESurfaceElement* pold = pme; 
					ss.m_data[i].m_rs = vec2d(0,0);

					pme = cpp.Project(m, q, ss.m_data[i].m_rs);

					if (pme == 0)
					{
						// nope, it has genuinly left contact
						int* n = &pold->m_node[0];
//						feLog("node %d has left element (%d, %d, %d, %d)\n", m+1, n[0]+1, n[1]+1, n[2]+1, n[3]+1);
					}
					else 
					{
/*						if (pme != pold)
						{
							feLog("node %d has switched segments: ", m + 1);
							int* n = &pold->m_node[0];
							feLog("from (%d, %d, %d, %d), ", n[0] + 1, n[1] + 1, n[2] + 1, n[3] + 1);
							n = &pme->m_node[0];
							feLog("to (%d, %d, %d, %d)\n", n[0] + 1, n[1] + 1, n[2] + 1, n[3] + 1);
						}
*/
						if (m_mu*m_epsf > 0)
						{
							// the node has moved to another segment.
							// If friction is active we need to translate the frictional
							// data to the new segment.
							FESurfaceElement& eo = *pold;
							FESurfaceElement& en = *pme;
							MapFrictionData(i, ss, ms, en, eo, q);
						}
					}
				}
			}
		}
		else if (bupseg)
		{
			// get the secondary surface element
			// don't forget to initialize the search for the first node!
			ss.m_data[i].m_rs = vec2d(0,0);
			pme = cpp.Project(m, q, ss.m_data[i].m_rs);
			if (pme)
			{
				// the node has come into contact so make sure to initialize
				// the previous natural coordinates for friction.
				ss.m_data[i].m_rsp = ss.m_data[i].m_rs;
			}
		}

		// if we found a secondary surface element, update the gap and normal data
		ss.m_data[i].m_pme = pme;
		if (pme != 0)
		{
			FESurfaceElement& mel =  *ss.m_data[i].m_pme;

			r = ss.m_data[i].m_rs[0];
			s = ss.m_data[i].m_rs[1];

			// if this is a new contact, copy the current coordinates
			// to the previous ones
			ss.m_data[i].m_M = ss.Metric0(mel, r, s);

			// the normal is set to the secondary surface element normal
			ss.m_data[i].m_nu = ss.SurfaceNormal(mel, r, s);

			// calculate gap
			ss.m_data[i].m_gap = -(ss.m_data[i].m_nu*(x - q)) + ss.m_data[i].m_off;
			if (bmove && (ss.m_data[i].m_gap>0))
			{
				node.m_r0 = node.m_rt = q + ss.m_data[i].m_nu*ss.m_data[i].m_off;
				ss.m_data[i].m_gap = 0;
			}

			// TODO: what should we do if the gap function becomes
			// negative? setting the Lagrange multipliers to zero
			// might make the system unstable.
/*			if (ss.gap[i] < 0)
			{
				ss.Lm[i] = 0;
				ss.Lt[i][0] = 0;
				ss.Lt[i][1] = 0;
				ss.pme[i] = 0;
			}
*/		}
		else
		{
			// TODO: Is this a good criteria for out-of-contact?
			//		 perhaps this is not even necessary.
			// since the node is not in contact, we set the gap function 
			// and Lagrangian multiplier to zero
			ss.m_data[i].m_gap = 0;
			ss.m_data[i].m_Lm  = 0;
			ss.m_data[i].m_Lt[0] = ss.m_data[i].m_Lt[1] = 0;
		}
	}
}

//-----------------------------------------------------------------------------
//! updates sliding interface data
//! niter is the number of Newton iterations.
void FESlidingInterface::Update()
{
	int niter = GetFEModel()->GetCurrentStep()->GetFESolver()->m_niter;

	// should we do a segment update or not?
	// TODO: check what happens when m_nsegup == -1 and m_npass = 2;
	// We have to make sure that in this case, both surfaces get at least
	// one pass!
	bool bupdate = (m_bfirst || (m_nsegup == 0)? true : (niter <= m_nsegup));

	// project primary surface onto secondary surface
	// this also calculates the nodal gap functions
	ProjectSurface(m_ss, m_ms, bupdate);
	if (m_btwo_pass) ProjectSurface(m_ms, m_ss, bupdate);

	// Update the net contact pressures
	UpdateContactPressures();

	// set the first-entry-flag to false
	m_bfirst = false;
}

//-----------------------------------------------------------------------------

void FESlidingInterface::LoadVector(FEGlobalVector& R, const FETimeInfo& tp)
{
	// element contact force vector
	vector<double> fe;

	// the lm array for this force vector
	vector<int> lm;

	// the en array
	vector<int> en;

	// the elements LM vectors
	vector<int> sLM;
	vector<int> mLM;

	const int MN = FEElement::MAX_NODES;
	vec3d r0[MN];
	double w[MN];
	double* Gr, *Gs;
	double detJ[MN];
	vec3d dxr, dxs;

	// do two-pass
	int npass = (m_btwo_pass?2:1);
	for (int np=0; np<npass; ++np)
	{
		// pick the primary and secondary surfaces
		FESlidingSurface& ss = (np==0? m_ss : m_ms);
		FESlidingSurface& ms = (np==0? m_ms : m_ss);

		// loop over all primary surface facets
		int ne = ss.Elements();
		for (int j=0; j<ne; ++j)
		{
			// get the next element
			FESurfaceElement& sel = ss.Element(j);
			int nseln = sel.Nodes();

			// get the element's LM array
			ss.UnpackLM(sel, sLM);

			// nodal coordinates
			for (int i=0; i<nseln; ++i) r0[i] = ss.GetMesh()->Node(sel.m_node[i]).m_r0;

			// we calculate all the metrics we need before we
			// calculate the nodal forces
			for (int n=0; n<nseln; ++n)
			{
				Gr = sel.Gr(n);
				Gs = sel.Gs(n);

				// calculate jacobian
				// note that we are integrating over the reference surface
				dxr = dxs = vec3d(0,0,0);
				for (int k=0; k<nseln; ++k)
				{
					dxr.x += Gr[k]*r0[k].x;
					dxr.y += Gr[k]*r0[k].y;
					dxr.z += Gr[k]*r0[k].z;

					dxs.x += Gs[k]*r0[k].x;
					dxs.y += Gs[k]*r0[k].y;
					dxs.z += Gs[k]*r0[k].z;
				}

				// jacobians
				detJ[n] = (dxr ^ dxs).norm();

				// integration weights
				w[n] = sel.GaussWeights()[n];
			}

			// loop over primary surface element nodes (which are the integration points as well)
			// and calculate the contact nodal force
			for (int n=0; n<nseln; ++n)
			{
				// get the local node number
				int m = sel.m_lnode[n];

				// see if this node's constraint is active
				// that is, if it has an element associated with it
				// TODO: is this a good way to test for an active constraint
				// The rigid wall criteria seems to work much better.
				if (ss.m_data[m].m_pme != 0)
				{
					// This node is active and could lead to a non-zero
					// contact force.
					// get the secondary surface element
					FESurfaceElement& mel = *ss.m_data[m].m_pme;
					ms.UnpackLM(mel, mLM);

					// calculate the degrees of freedom
					int nmeln = mel.Nodes();
					int ndof = 3*(nmeln+1);
					fe.resize(ndof);

					// calculate the nodal force
					ContactNodalForce(m, ss, mel, fe);

					// multiply force with weights
					for (int l=0; l<ndof; ++l) fe[l] *= detJ[n]*w[n];
					
					// fill the lm array
					lm.resize(3*(nmeln+1));
					lm[0] = sLM[n*3  ];
					lm[1] = sLM[n*3+1];
					lm[2] = sLM[n*3+2];

					for (int l=0; l<nmeln; ++l)
					{
						lm[3*(l+1)  ] = mLM[l*3  ];
						lm[3*(l+1)+1] = mLM[l*3+1];
						lm[3*(l+1)+2] = mLM[l*3+2];
					}

					// fill the en array
					en.resize(nmeln+1);
					en[0] = sel.m_node[n];
					for (int l=0; l<nmeln; ++l) en[l+1] = mel.m_node[l];

					// assemble into global force vector
					R.Assemble(en, lm, fe);
				}
			}
		}
	}
}

//-----------------------------------------------------------------------------
//! Calculates the contact force on a node.
//! \param[in] m local node number
//! \param[out] fe force vector

void FESlidingInterface::ContactNodalForce(int m, FESlidingSurface& ss, FESurfaceElement& mel, vector<double>& fe)
{
	vec3d dxr, dxs;

	// normal force
	double tn, Ln;

	// gap function
	double gap;

	// tangents
	vec3d tau1, tau2;

	// max nr of element nodes
	const int MAXMN = FEElement::MAX_NODES;

	// secondary surface element nodes
	vec3d rtm[MAXMN];

	// shape function values
	double H[MAXMN], Hr[MAXMN], Hs[MAXMN];

	// contact vectors
	double N[3*(MAXMN+1)], N1[3*(MAXMN+1)], N2[3*(MAXMN+1)];
	double T1[3*(MAXMN+1)], T2[3*(MAXMN+1)], D1[3*(MAXMN+1)], D2[3*(MAXMN+1)];

	// surface metrics
	double A[2][2], M[2][2], K[2][2];
	double detA;

	double eps, scale = Penalty();

	// get the mesh
	FEMesh& mesh = GetFEModel()->GetMesh();

	double Tt[2];

	int nmeln, ndof;

	// gap function
	gap = ss.m_data[m].m_gap;

	// normal penalty
	eps = ss.m_data[m].m_eps*scale;

	// get primary surface node normal force
	Ln = ss.m_data[m].m_Lm;
	tn = Ln + eps*gap;
	tn = MBRACKET(tn);

	// get the primary surface  node normal
	vec3d& nu = ss.m_data[m].m_nu;

	nmeln = mel.Nodes();
	ndof = 3*(1 + nmeln);

	// metric tensors
	mat2d Mk = ss.m_data[m].m_M;
	mat2d Mki = Mk.inverse();

	// get the secondary surface element node positions
	for (int k=0; k<nmeln; ++k) rtm[k] = mesh.Node(mel.m_node[k]).m_rt;

	// isoparametric coordinates of the projected node
	// onto the secondary surface element
	double r = ss.m_data[m].m_rs[0];
	double s = ss.m_data[m].m_rs[1];

	// get the coordinates at the previous step
	double rp = ss.m_data[m].m_rsp[0];
	double sp = ss.m_data[m].m_rsp[1];

	// get the secondary surface element shape function values at this node
	mel.shape_fnc(H, r, s);

	// --- N O R M A L   T R A C T I O N ---

	// calculate contact vectors for normal traction
	N[0] = nu.x;
	N[1] = nu.y;
	N[2] = nu.z;
	for (int l=0; l<nmeln; ++l)
	{
		N[3*(l+1)  ] = -H[l]*nu.x;
		N[3*(l+1)+1] = -H[l]*nu.y;
		N[3*(l+1)+2] = -H[l]*nu.z;
	}

	// calculate force vector
	for (int l=0; l<ndof; ++l) fe[l] = tn*N[l];

	// --- T A N G E N T I A L   T R A C T I O N ---
	if (m_mu*m_epsf > 0)
	{
		// Lagrangian traction
		double Lt[2];
		Lt[0] = ss.m_data[m].m_Lt[0];
		Lt[1] = ss.m_data[m].m_Lt[1];

		// calculate contact vector for tangential traction
		// only if both the friction coefficient and friction
		// penalty factor are non-zero

		// get the secondary surface shape function derivative values at this node
		mel.shape_deriv(Hr, Hs, r, s);

		// get the tangent vectors
		tau1 = tau2 = vec3d(0,0,0);
		for (int k=0; k<nmeln; ++k)
		{
			tau1.x += Hr[k]*rtm[k].x;
			tau1.y += Hr[k]*rtm[k].y;
			tau1.z += Hr[k]*rtm[k].z;
		
			tau2.x += Hs[k]*rtm[k].x;
			tau2.y += Hs[k]*rtm[k].y;
			tau2.z += Hs[k]*rtm[k].z;
		}

		// set up the Ti vectors
		T1[0] = tau1.x; T2[0] = tau2.x;
		T1[1] = tau1.y; T2[1] = tau2.y;
		T1[2] = tau1.z; T2[2] = tau2.z;

		for (int k=0; k<nmeln; ++k)
		{
			T1[(k+1)*3  ] = -H[k]*tau1.x;
			T1[(k+1)*3+1] = -H[k]*tau1.y;
			T1[(k+1)*3+2] = -H[k]*tau1.z;

			T2[(k+1)*3  ] = -H[k]*tau2.x;
			T2[(k+1)*3+1] = -H[k]*tau2.y;
			T2[(k+1)*3+2] = -H[k]*tau2.z;
		}

		// set up the Ni vectors
		N1[0] = N2[0] = 0;
		N1[1] = N2[1] = 0;
		N1[2] = N2[2] = 0;

		for (int k=0; k<nmeln; ++k)
		{
			N1[(k+1)*3  ] = -Hr[k]*nu.x;
			N1[(k+1)*3+1] = -Hr[k]*nu.y;
			N1[(k+1)*3+2] = -Hr[k]*nu.z;

			N2[(k+1)*3  ] = -Hs[k]*nu.x;
			N2[(k+1)*3+1] = -Hs[k]*nu.y;
			N2[(k+1)*3+2] = -Hs[k]*nu.z;
		}

		// calculate metric tensor
		M[0][0] = tau1*tau1; M[0][1] = tau1*tau2; 
		M[1][0] = tau2*tau1; M[1][1] = tau2*tau2; 

		// calculate curvature tensor
		K[0][0] = 0; K[0][1] = 0;
		K[1][0] = 0; K[1][1] = 0;

		double Grr[FEElement::MAX_NODES];
		double Grs[FEElement::MAX_NODES];
		double Gss[FEElement::MAX_NODES];
		mel.shape_deriv2(Grr, Grs, Gss, r, s);
		for (int k=0; k<nmeln; ++k)
		{
			K[0][0] += (nu*rtm[k])*Grr[k];
			K[0][1] += (nu*rtm[k])*Grs[k];
			K[1][0] += (nu*rtm[k])*Grs[k];
			K[1][1] += (nu*rtm[k])*Gss[k];
		}

		// setup A matrix
		A[0][0] = M[0][0] + gap*K[0][0];
		A[0][1] = M[0][1] + gap*K[0][1];
		A[1][0] = M[1][0] + gap*K[1][0];
		A[1][1] = M[1][1] + gap*K[1][1];

		detA = A[0][0]*A[1][1] - A[0][1]*A[1][0];

		// setup Di vectors
		for (int k=0; k<ndof; ++k)
		{
			D1[k] = (1/detA)*(A[1][1]*(T1[k]+gap*N1[k]) - A[0][1]*(T2[k] + gap*N2[k]));
			D2[k] = (1/detA)*(A[0][0]*(T2[k]+gap*N2[k]) - A[0][1]*(T1[k] + gap*N1[k]));
		}

		// calculate friction tractions
		// a. calculate trial state
		Tt[0] = Lt[0] + m_epsf*(Mk[0][0]*(r - rp) + Mk[0][1]*(s - sp));
		Tt[1] = Lt[1] + m_epsf*(Mk[1][0]*(r - rp) + Mk[1][1]*(s - sp));

		double TMT = Tt[0]*(Mki[0][0]*Tt[0]+Mki[0][1]*Tt[1])+Tt[1]*(Mki[1][0]*Tt[0]+Mki[1][1]*Tt[1]);
		assert(TMT >= 0);

		double phi = sqrt(TMT) - m_mu* tn;

		// b. return map
		if (phi > 0)
		{
			Tt[0] = m_mu* tn *Tt[0]/sqrt(TMT);
			Tt[1] = m_mu* tn *Tt[1]/sqrt(TMT);
		}

		// tangential force vector
		for (int l=0; l<ndof; ++l) fe[l] -= (Tt[0]*D1[l] + Tt[1]*D2[l]);
	}
}

//-----------------------------------------------------------------------------

void FESlidingInterface::StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp)
{
	FEElementMatrix ke;

	const int MAXMN = FEElement::MAX_NODES;
	vector<int> lm(3*(MAXMN + 1));
	vector<int> en(MAXMN+1);

	double *Gr, *Gs, w[MAXMN];
	vec3d r0[MAXMN];

	double detJ[MAXMN];
	vec3d dxr, dxs;

	vector<int> sLM;
	vector<int> mLM;

	// do two-pass
	int npass = (m_btwo_pass?2:1);
	for (int np=0; np<npass; ++np)
	{
		// get the primary and secondary surface
		FESlidingSurface& ss = (np==0?m_ss:m_ms);	
		FESlidingSurface& ms = (np==0?m_ms:m_ss);	

		// loop over all primary surface elements
		int ne = ss.Elements();
		for (int j=0; j<ne; ++j)
		{
			// unpack the next element
			FESurfaceElement& se = ss.Element(j);
			int nseln = se.Nodes();

			// get the element's LM array
			ss.UnpackLM(se, sLM);

			// get the nodal coordinates
			for (int i=0; i<nseln; ++i) r0[i] = ss.GetMesh()->Node(se.m_node[i]).m_r0;

			// get all the metrics we need 
			for (int n=0; n<nseln; ++n)
			{
				Gr = se.Gr(n);
				Gs = se.Gs(n);

				// calculate jacobian
				dxr = dxs = vec3d(0,0,0);
				for (int k=0; k<nseln; ++k)
				{
					dxr.x += Gr[k]*r0[k].x;
					dxr.y += Gr[k]*r0[k].y;
					dxr.z += Gr[k]*r0[k].z;

					dxs.x += Gs[k]*r0[k].x;
					dxs.y += Gs[k]*r0[k].y;
					dxs.z += Gs[k]*r0[k].z;
				}

				detJ[n] = (dxr ^ dxs).norm();
				w[n] = se.GaussWeights()[n];
			}

			// loop over all integration points (that is nodes)
			for (int n=0; n<nseln; ++n)
			{
				int m = se.m_lnode[n];

				// see if this node's constraint is active
				// that is, if it has an element associated with it
				if (ss.m_data[m].m_pme != 0)
				{
					// get the secondary surface element
					FESurfaceElement& me = *ss.m_data[m].m_pme;

					// get the secondary surface element's LM array
					ms.UnpackLM(me, mLM);

					int nmeln = me.Nodes();
					int ndof = 3*(nmeln+1);

					// calculate the stiffness matrix
					ke.resize(ndof, ndof);
					ContactNodalStiffness(m, ss, me, ke);

					// muliply with weights
					for (int k=0; k<ndof; ++k)
						for (int l=0; l<ndof; ++l) ke[k][l] *= detJ[n]*w[n];

					// fill the lm array
					lm[0] = sLM[n*3  ];
					lm[1] = sLM[n*3+1];
					lm[2] = sLM[n*3+2];

					for (int k=0; k<nmeln; ++k)
					{
						lm[3*(k+1)  ] = mLM[k*3  ];
						lm[3*(k+1)+1] = mLM[k*3+1];
						lm[3*(k+1)+2] = mLM[k*3+2];
					}

					// create the en array
					en.resize(nmeln+1);
					en[0] = se.m_node[n];
					for (int k=0; k<nmeln; ++k) en[k+1] = me.m_node[k];
						
					// assemble stiffness matrix
					ke.SetNodes(en);
					ke.SetIndices(lm);
					LS.Assemble(ke);
				}
			}
		}
	}
}

//-----------------------------------------------------------------------------

void FESlidingInterface::ContactNodalStiffness(int m, FESlidingSurface& ss, FESurfaceElement& mel, matrix& ke)
{
	const int MAXMN = FEElement::MAX_NODES;

	vector<int> lm(3*(MAXMN+1));
	vector<int> en(MAXMN + 1);

	vec3d dxr, dxs;
	double H[MAXMN], Hr[MAXMN], Hs[MAXMN];

	double N[3*(MAXMN+1)], T1[3*(MAXMN+1)], T2[3*(MAXMN+1)];
	double N1[3*(MAXMN+1)], N2[3*(MAXMN+1)], D1[3*(MAXMN+1)], D2[3*(MAXMN+1)];
	double Nb1[3*(MAXMN+1)], Nb2[3*(MAXMN+1)];

	// get the mesh
	FEMesh& mesh = GetFEModel()->GetMesh();

	// nr of element nodes and degrees of freedom 
	int nmeln = mel.Nodes();
	int ndof = 3*(1 + nmeln);

	// penalty factor
	double scale = Penalty();
	double eps = ss.m_data[m].m_eps*scale;

	// nodal coordinates
	vec3d rt[MAXMN];
	for (int j=0; j<nmeln; ++j) rt[j] = mesh.Node(mel.m_node[j]).m_rt;

	// node natural coordinates in secondary surface element
	double r = ss.m_data[m].m_rs[0];
	double s = ss.m_data[m].m_rs[1];

	// gap
	double gap = ss.m_data[m].m_gap;

	// lagrange multiplier
	double Lm = ss.m_data[m].m_Lm;

	// get node normal force
	double tn = Lm + eps*gap;
	tn = MBRACKET(tn);

	// get the node normal
	vec3d& nu = ss.m_data[m].m_nu;

	// get the secondary surface element shape function values and the derivatives at this node
	mel.shape_fnc(H, r, s);
	mel.shape_deriv(Hr, Hs, r, s);

	// get the tangent vectors
	vec3d tau[2];
	ss.CoBaseVectors(mel, r, s, tau);

	// set up the N vector
	N[0] = nu.x;
	N[1] = nu.y;
	N[2] = nu.z;

	for (int k=0; k<nmeln; ++k)
	{
		N[(k+1)*3  ] = -H[k]*nu.x;
		N[(k+1)*3+1] = -H[k]*nu.y;
		N[(k+1)*3+2] = -H[k]*nu.z;
	}
	
	// set up the Ti vectors
	T1[0] = tau[0].x; T2[0] = tau[1].x;
	T1[1] = tau[0].y; T2[1] = tau[1].y;
	T1[2] = tau[0].z; T2[2] = tau[1].z;

	for (int k=0; k<nmeln; ++k)
	{
		T1[(k+1)*3  ] = -H[k]*tau[0].x;
		T1[(k+1)*3+1] = -H[k]*tau[0].y;
		T1[(k+1)*3+2] = -H[k]*tau[0].z;

		T2[(k+1)*3  ] = -H[k]*tau[1].x;
		T2[(k+1)*3+1] = -H[k]*tau[1].y;
		T2[(k+1)*3+2] = -H[k]*tau[1].z;
	}

	// set up the Ni vectors
	N1[0] = N2[0] = 0;
	N1[1] = N2[1] = 0;
	N1[2] = N2[2] = 0;

	for (int k=0; k<nmeln; ++k)
	{
		N1[(k+1)*3  ] = -Hr[k]*nu.x;
		N1[(k+1)*3+1] = -Hr[k]*nu.y;
		N1[(k+1)*3+2] = -Hr[k]*nu.z;

		N2[(k+1)*3  ] = -Hs[k]*nu.x;
		N2[(k+1)*3+1] = -Hs[k]*nu.y;
		N2[(k+1)*3+2] = -Hs[k]*nu.z;
	}

	// calculate metric tensor
	mat2d M;
	M[0][0] = tau[0]*tau[0]; M[0][1] = tau[0]*tau[1]; 
	M[1][0] = tau[1]*tau[0]; M[1][1] = tau[1]*tau[1]; 

	// calculate reciprocal metric tensor
	mat2d Mi = M.inverse();

	// calculate curvature tensor
	double K[2][2] = {0};
	double Grr[FEElement::MAX_NODES];
	double Grs[FEElement::MAX_NODES];
	double Gss[FEElement::MAX_NODES];
	mel.shape_deriv2(Grr, Grs, Gss, r, s);
	for (int k=0; k<nmeln; ++k)
	{
		K[0][0] += (nu*rt[k])*Grr[k];
		K[0][1] += (nu*rt[k])*Grs[k];
		K[1][0] += (nu*rt[k])*Grs[k];
		K[1][1] += (nu*rt[k])*Gss[k];
	}

	// setup A matrix A = M + gK
	double A[2][2];
	A[0][0] = M[0][0] + gap*K[0][0];
	A[0][1] = M[0][1] + gap*K[0][1];
	A[1][0] = M[1][0] + gap*K[1][0];
	A[1][1] = M[1][1] + gap*K[1][1];

	// calculate determinant of A
	double detA = A[0][0]*A[1][1] - A[0][1]*A[1][0];

	// setup Di vectors
	for (int k=0; k<ndof; ++k)
	{
		D1[k] = (1/detA)*(A[1][1]*(T1[k]+gap*N1[k]) - A[0][1]*(T2[k] + gap*N2[k]));
		D2[k] = (1/detA)*(A[0][0]*(T2[k]+gap*N2[k]) - A[0][1]*(T1[k] + gap*N1[k]));
	}

	// setup Nbi vectors
	for (int k=0; k<ndof; ++k)
	{
		Nb1[k] = N1[k] - K[0][1]*D2[k];
		Nb2[k] = N2[k] - K[0][1]*D1[k];
	}

	// --- N O R M A L   S T I F F N E S S ---
	double sum;
	for (int k=0; k<ndof; ++k)
		for (int l=0; l<ndof; ++l)
			{
				sum = 0;

				sum = Mi[0][0]*Nb1[k]*Nb1[l]+Mi[0][1]*(Nb1[k]*Nb2[l]+Nb2[k]*Nb1[l])+Mi[1][1]*Nb2[k]*Nb2[l];
				sum *= gap;
				sum -= D1[k]*N1[l]+D2[k]*N2[l]+N1[k]*D1[l]+N2[k]*D2[l];
				sum += K[0][1]*(D1[k]*D2[l]+D2[k]*D1[l]);
				sum *= tn*m_knmult;

				sum += eps*HEAVYSIDE(Lm+eps*gap)*N[k]*N[l];
	
				ke[k][l] = sum;
			}

	// --- T A N G E N T I A L   S T I F F N E S S ---
	// We only calculate the tangential stiffness if friction is enabled. We also
	// make sure that the gap >= 0, i.e. the node is actually in contact, otherwise
	// I've noticed that the solution can diverge quickly.
	if ((m_mu*m_epsf > 0) && (gap >=0))
	{
		// get the traction multipliers
		double Lt[2];
		Lt[0] = ss.m_data[m].m_Lt[0];
		Lt[1] = ss.m_data[m].m_Lt[1];

		// get the metric tensor and its inverse
		mat2d& Mk = ss.m_data[m].m_M;
		mat2d Mki = Mk.inverse();

		// get the previous isoparameteric coordinates
		double rp = ss.m_data[m].m_rsp[0];
		double sp = ss.m_data[m].m_rsp[1];

		// get the traction
		double Tt[2];
		// a. trial state
		Tt[0] = Lt[0] + Mk[0][0]*(r - rp) + Mk[0][1]*(s - sp);
		Tt[1] = Lt[1] + Mk[1][0]*(r - rp) + Mk[1][1]*(s - sp);

		double TMT = Tt[0]*(Mki[0][0]*Tt[0] + Mki[0][1]*Tt[1]) + Tt[1]*(Mki[1][0]*Tt[0] + Mki[1][1]*Tt[1]);

		// calculate the normalized traction vector
		vec3d pt = tau[0]*Tt[0] + tau[1]*Tt[1];
		pt.unit();
		
		// calculate the covariant version
		double Pt[2] = {Tt[0]/sqrt(TMT), Tt[1]/sqrt(TMT)};
		double Ptc[2];
		Ptc[0] = Mki[0][0]*Pt[0]+Mki[0][1]*Pt[1];
		Ptc[1] = Mki[1][0]*Pt[0]+Mki[1][1]*Pt[1];

		//b. return map
		bool bstick = true;
		double phi = sqrt(TMT) - m_mu*tn;
		if (phi > 0)
		{
			Tt[0] = m_mu*Tt[0]/sqrt(TMT);
			Tt[1] = m_mu*Tt[1]/sqrt(TMT);
			bstick = false;
		}

		// we need to define additional arrays for the tangential
		// contribution of the contact stiffness
		double T11[3*(MAXMN+1)] = {0}, T12[3*(MAXMN+1)] = {0}, T21[3*(MAXMN+1)] = {0}, T22[3*(MAXMN+1)] = {0};	// Tab matrix
		double N11[3*(MAXMN+1)] = {0}, N12[3*(MAXMN+1)] = {0}, N21[3*(MAXMN+1)] = {0}, N22[3*(MAXMN+1)] = {0};	// Nab matrix
		double P1[3*(MAXMN+1)] = {0}, P2[3*(MAXMN+1)] = {0};	// P arrays
		double Tb11[3*(MAXMN+1)], Tb21[3*(MAXMN+1)], Tb12[3*(MAXMN+1)], Tb22[3*(MAXMN+1)]; // Tbar matrix
		double Pb1[3*(MAXMN+1)], Pb2[3*(MAXMN+1)]; // Pbar array

		for (int k=0; k<nmeln; ++k)
		{
			T11[3*(k+1)  ] = -Hr[k]*tau[0].x;
			T11[3*(k+1)+1] = -Hr[k]*tau[0].y;
			T11[3*(k+1)+2] = -Hr[k]*tau[0].z;

			T12[3*(k+1)  ] = -Hs[k]*tau[0].x;
			T12[3*(k+1)+1] = -Hs[k]*tau[0].y;
			T12[3*(k+1)+2] = -Hs[k]*tau[0].z;

			T21[3*(k+1)  ] = -Hr[k]*tau[1].x;
			T21[3*(k+1)+1] = -Hr[k]*tau[1].y;
			T21[3*(k+1)+2] = -Hr[k]*tau[1].z;

			T22[3*(k+1)  ] = -Hs[k]*tau[1].x;
			T22[3*(k+1)+1] = -Hs[k]*tau[1].y;
			T22[3*(k+1)+2] = -Hs[k]*tau[1].z;

			if (nmeln==4)
			{
				N12[3*(k+1)  ] = N21[3*(k+1)  ] = -0.25*nu.x;
				N12[3*(k+1)+1] = N21[3*(k+1)+1] = -0.25*nu.y;
				N12[3*(k+1)+2] = N21[3*(k+1)+2] = -0.25*nu.z;
			}
			else if (nmeln == 6) assert(false);

			P1[3*(k+1)  ] = -Hr[k]*pt.x;
			P1[3*(k+1)+1] = -Hr[k]*pt.y;
			P1[3*(k+1)+2] = -Hr[k]*pt.z;

			P2[3*(k+1)  ] = -Hs[k]*pt.x;
			P2[3*(k+1)+1] = -Hs[k]*pt.y;
			P2[3*(k+1)+2] = -Hs[k]*pt.z;
		}

		vec3d g12(0,0,0);
		if (nmeln == 4) 
		{
			const double Grs[4] = {0.25, -0.25, 0.25, -0.25};
			g12 = rt[0]*Grs[0]+rt[1]*Grs[1]+rt[2]*Grs[2]+rt[3]*Grs[3];
		}
		else if (nmeln == 6)
		{
			const double Grs[6] = {4.0, 0.0, 0.0, -4.0, 4.0, -4.0};
			g12 = rt[0]*Grs[0]+rt[1]*Grs[1]+rt[2]*Grs[2]+rt[3]*Grs[3]+rt[4]*Grs[4]+rt[5]*Grs[5];
		}
		else
		{
			assert(false);
		}
		double gt1 = g12*tau[0];
		double gt2 = g12*tau[1];
		double gp = g12*pt;

		for (int k=0; k<ndof; ++k)
		{
			Tb11[k] = T11[k] - gt1*D2[k];
			Tb12[k] = T12[k] - gt1*D1[k];

			Tb21[k] = T21[k] - gt2*D2[k];
			Tb22[k] = T22[k] - gt2*D1[k];

			Pb1[k] = P1[k] - gp*D2[k];
			Pb2[k] = P2[k] - gp*D1[k];
		}

		// raise the indices of A
		double Ac[2][2];
		for (int k=0; k<2; ++k)
			for (int l=0; l<2; ++l)
			{
				Ac[k][l] = 0;
				for (int i=0; i<2; ++i)
					for (int j=0; j<2; ++j) Ac[k][l] += Mki[k][i]*Mki[l][j]*A[i][j];
			}

		vec3d Hrs[2][2] = {{vec3d(0,0,0),vec3d(0,0,0)},{vec3d(0,0,0),vec3d(0,0,0)}};
		if (nmeln == 4) 
		{
			const double Grs[4] = {0.25, -0.25, 0.25, -0.25};
			Hrs[0][1] = Hrs[1][0] = rt[0]*Grs[0]+rt[1]*Grs[1]+rt[2]*Grs[2]+rt[3]*Grs[3];
		}
		else if (nmeln == 6)
		{
			const double Grs[6] = {4.0, 0.0, 0.0, -4.0, 4.0, -4.0};
			Hrs[0][1] = Hrs[1][0] = rt[0]*Grs[0]+rt[1]*Grs[1]+rt[2]*Grs[2]+rt[3]*Grs[3]+rt[4]*Grs[4]+rt[5]*Grs[5];
		}

		// calculate stiffness matrix
		// NOTE: I think I need Mi (iso Mki) for KT1 and KT2. I only need Mk (iso M) for the direct stiffnesses
		double kij;
		for (int i=0; i<ndof; ++i)
			for (int j=0; j<ndof; ++j)
			{
				// KT1
				kij  = T11[i]*D1[j] + T12[i]*D2[j];
				kij += D1[i]*T11[j] + D2[i]*T12[j];
				kij -= (Hrs[0][1]*tau[0])*D1[i]*D2[j] + (Hrs[1][0]*tau[0])*D2[i]*D1[j];
				kij += Tb11[i]*D1[j] + Tb21[i]*D2[j];
				kij += D1[i]*Tb11[j] + D2[i]*Tb21[j];
				kij += gap*(N11[i]*D1[j] + N12[i]*D2[j] + D1[i]*N11[j] + D2[i]*N12[j]);
				kij -= N[i]*Nb1[j] - Nb1[i]*N[j];
				kij -= T1[i]*Mi[0][0]*Tb11[j] + T1[i]*Mi[0][1]*Tb21[j] + T2[i]*Mi[1][0]*Tb11[j] + T2[i]*Mi[1][1]*Tb21[j];
				kij -= Tb11[i]*Mi[0][0]*T1[j] + Tb21[i]*Mi[0][1]*T1[j] + Tb11[i]*Mi[1][0]*T2[j] + Tb21[i]*Mi[1][1]*T2[j];

				ke[i][j] += m_ktmult*(Tt[0]*Ac[0][0] + Tt[1]*Ac[1][0])*kij;

				// KT2
				kij  = T21[i]*D1[j] + T22[i]*D2[j];
				kij += D1[i]*T21[j] + D2[i]*T22[j];
				kij -= (Hrs[0][1]*tau[1])*D1[i]*D2[j] + (Hrs[1][0]*tau[1])*D2[i]*D1[j];
				kij += Tb12[i]*D1[j] + Tb22[i]*D2[j];
				kij += D1[i]*Tb12[j] + D2[i]*Tb22[j];
				kij += gap*(N21[i]*D1[j] + N22[i]*D2[j] + D1[i]*N21[j] + D2[i]*N22[j]);
				kij -= N[i]*Nb2[j] - Nb2[i]*N[j];
				kij -= T1[i]*Mi[0][0]*Tb12[j] + T1[i]*Mi[0][1]*Tb22[j] + T2[i]*Mi[1][0]*Tb12[j] + T2[i]*Mi[1][1]*Tb22[j];
				kij -= Tb12[i]*Mi[0][0]*T1[j] + Tb22[i]*Mi[0][1]*T1[j] + Tb12[i]*Mi[1][0]*T2[j] + Tb22[i]*Mi[1][1]*T2[j];

				ke[i][j] += m_ktmult*(Tt[0]*Ac[0][1] + Tt[1]*Ac[1][1])*kij;

				// kdirect
				if (bstick)
				{
					kij = Mk[0][0]*D1[i]*D1[j] + Mk[0][1]*D1[i]*D2[j] + Mk[1][0]*D2[i]*D1[j] + Mk[1][1]*D2[i]*D2[j];
					ke[i][j] += m_epsf*kij;
				}
				else
				{
					kij  = (1.0 - Ptc[0]*Pt[0])*(Mk[0][0]*D1[i]*D1[j]+Mk[0][1]*D1[i]*D2[j]);
					kij += (    - Ptc[0]*Pt[1])*(Mk[1][0]*D1[i]*D1[j]+Mk[1][1]*D1[i]*D2[j]);
					kij += (    - Ptc[1]*Pt[0])*(Mk[0][0]*D2[i]*D1[j]+Mk[0][1]*D2[i]*D2[j]);
					kij += (1.0 - Ptc[1]*Pt[1])*(Mk[1][0]*D2[i]*D1[j]+Mk[1][1]*D2[i]*D2[j]);
					
					ke[i][j] += m_ktmult*m_epsf*m_mu*tn/sqrt(TMT)*kij;
				}
			}
	}
}

//-----------------------------------------------------------------------------

bool FESlidingInterface::Augment(int naug, const FETimeInfo& tp)
{
	// make sure we need to augment
	if (m_laugon != 1) return true;

	double Ln;
	double Lt[2];
	bool bconv = true;
	mat2d Mi;

	// penalty factor
	double eps, scale = Penalty();

	// --- c a l c u l a t e   i n i t i a l   n o r m s ---
	// a. normal component
	double normL0 = 0;
	for (int i=0; i<m_ss.Nodes(); ++i)	normL0 += m_ss.m_data[i].m_Lm*m_ss.m_data[i].m_Lm;
	for (int i=0; i<m_ms.Nodes(); ++i)	normL0 += m_ms.m_data[i].m_Lm*m_ms.m_data[i].m_Lm;

	// b. tangential component
	if (m_mu*m_epsf > 0)
	{
		for (int i=0; i<m_ss.Nodes(); ++i)
		{
			if (m_ss.m_data[i].m_pme)
			{
				Lt[0] = m_ss.m_data[i].m_Lt[0];
				Lt[1] = m_ss.m_data[i].m_Lt[1];
				mat2d& M = m_ss.m_data[i].m_M;
				Mi = M.inverse();
				normL0 += Lt[0]*(Mi[0][0]*Lt[0] + Mi[0][1]*Lt[1]) + Lt[1]*(Mi[1][0]*Lt[0] + Mi[1][1]*Lt[1]);
			}
		}

		for (int i=0; i<m_ms.Nodes(); ++i)
		{
			if (m_ms.m_data[i].m_pme)
			{
				Lt[0] = m_ms.m_data[i].m_Lt[0];
				Lt[1] = m_ms.m_data[i].m_Lt[1];
				mat2d& M = m_ms.m_data[i].m_M;
				Mi = M.inverse();
				normL0 += Lt[0]*(Mi[0][0]*Lt[0] + Mi[0][1]*Lt[1]) + Lt[1]*(Mi[1][0]*Lt[0] + Mi[1][1]*Lt[1]);
			}
		}
	}
	normL0 = sqrt(normL0);

	// --- c a l c u l a t e   c u r r e n t   n o r m s ---
	// a. normal component
	double normL1 = 0;	// force norm
	double normg1 = 0;	// gap norm
	int N = 0;
	for (int i=0; i<m_ss.Nodes(); ++i)
	{
		eps = m_ss.m_data[i].m_eps*scale;

		// update Lagrange multipliers
		Ln = m_ss.m_data[i].m_Lm + eps*m_ss.m_data[i].m_gap;
		Ln = MBRACKET(Ln);

		normL1 += Ln*Ln;

		if (m_ss.m_data[i].m_gap > 0)
		{
			normg1 += m_ss.m_data[i].m_gap*m_ss.m_data[i].m_gap;
			++N;
		}
	}	

	for (int i=0; i<m_ms.Nodes(); ++i)
	{
		eps = m_ms.m_data[i].m_eps*scale;

		// update Lagrange multipliers
		Ln = m_ms.m_data[i].m_Lm + eps*m_ms.m_data[i].m_gap;
		Ln = MBRACKET(Ln);

		normL1 += Ln*Ln;
		if (m_ms.m_data[i].m_gap > 0)
		{
			normg1 += m_ms.m_data[i].m_gap*m_ms.m_data[i].m_gap;
			++N;
		}
	}
	if (N == 0) N=1;

	// b. tangential component
	if (m_mu*m_epsf > 0)
	{
		double r, s, rp, sp;
		for (int i=0; i<m_ss.Nodes(); ++i)
		{
			if (m_ss.m_data[i].m_pme)
			{
				r = m_ss.m_data[i].m_rs[0];
				s = m_ss.m_data[i].m_rs[1];
				rp = m_ss.m_data[i].m_rsp[0];
				sp = m_ss.m_data[i].m_rsp[1];
				Ln = m_ss.m_data[i].m_Lm;
	
				mat2d& Mk = m_ss.m_data[i].m_M;
				Mi = Mk.inverse();

				Lt[0] = m_ss.m_data[i].m_Lt[0] + m_epsf*(Mk[0][0]*(r - rp) + Mk[0][1]*(s - sp));
				Lt[1] = m_ss.m_data[i].m_Lt[1] + m_epsf*(Mk[1][0]*(r - rp) + Mk[1][1]*(s - sp));

				double TMT = Lt[0]*(Mi[0][0]*Lt[0]+Mi[0][1]*Lt[1])+Lt[1]*(Mi[1][0]*Lt[0]+Mi[1][1]*Lt[1]);
				double phi = sqrt(TMT) - m_mu*Ln;

				// b. return map
				if (phi > 0)
				{
					Lt[0] = m_mu*Ln*Lt[0]/sqrt(TMT);
					Lt[1] = m_mu*Ln*Lt[1]/sqrt(TMT);
				}

				normL1 += Lt[0]*(Mi[0][0]*Lt[0] + Mi[0][1]*Lt[1]) + Lt[1]*(Mi[1][0]*Lt[0] + Mi[1][1]*Lt[1]);
			}
		}

		for (int i=0; i<m_ms.Nodes(); ++i)
		{
			if (m_ms.m_data[i].m_pme)
			{
				r = m_ms.m_data[i].m_rs[0];
				s = m_ms.m_data[i].m_rs[1];
				rp = m_ms.m_data[i].m_rsp[0];
				sp = m_ms.m_data[i].m_rsp[1];
				Ln = m_ms.m_data[i].m_Lm;

				mat2d& Mk = m_ms.m_data[i].m_M;
				Mi = Mk.inverse();

				Lt[0] = m_ms.m_data[i].m_Lt[0] + m_epsf*(Mk[0][0]*(r - rp) + Mk[0][1]*(s - sp));
				Lt[1] = m_ms.m_data[i].m_Lt[1] + m_epsf*(Mk[1][0]*(r - rp) + Mk[1][1]*(s - sp));

				double TMT = Lt[0]*(Mi[0][0]*Lt[0]+Mi[0][1]*Lt[1])+Lt[1]*(Mi[1][0]*Lt[0]+Mi[1][1]*Lt[1]);
				double phi = sqrt(TMT) - m_mu*Ln;

				// b. return map
				if (phi > 0)
				{
					Lt[0] = m_mu*Ln*Lt[0]/sqrt(TMT);
					Lt[1] = m_mu*Ln*Lt[1]/sqrt(TMT);
				}

				normL1 += Lt[0]*(Mi[0][0]*Lt[0] + Mi[0][1]*Lt[1]) + Lt[1]*(Mi[1][0]*Lt[0] + Mi[1][1]*Lt[1]);
			}
		}
	}

	normL1 = sqrt(normL1);
	normg1 = sqrt(normg1 / N);

	if (naug == 0) m_normg0 = 0;

	// calculate and print convergence norms
	double lnorm = 0, gnorm = 0;
	if (normL1 != 0) lnorm = fabs(normL1 - normL0)/normL1; else lnorm = fabs(normL1 - normL0);
	if (normg1 != 0) gnorm = fabs(normg1 - m_normg0)/normg1; else gnorm = fabs(normg1 - m_normg0);

	feLog(" sliding interface # %d\n", GetID());
	feLog("                        CURRENT        REQUIRED\n");
	feLog("    normal force : %15le", lnorm);
	if (m_atol > 0) feLog("%15le\n", m_atol); else feLog("       ***\n");
	feLog("    gap function : %15le", gnorm);
	if (m_gtol > 0) feLog("%15le\n", m_gtol); else feLog("       ***\n");

	// check convergence
	bconv = true;
	if ((m_atol > 0) && (lnorm > m_atol)) bconv = false;
	if ((m_gtol > 0) && (gnorm > m_gtol)) bconv = false;
	if (m_naugmin > naug) bconv = false;
	if (m_naugmax <= naug) bconv = true;
		
	if (bconv == false)
	{
		// we did not converge so update multipliers
		for (int i=0; i<m_ss.Nodes(); ++i)
		{
			eps = m_ss.m_data[i].m_eps*scale;

			// update Lagrange multipliers
			Ln = m_ss.m_data[i].m_Lm + eps*m_ss.m_data[i].m_gap;
			m_ss.m_data[i].m_Lm = MBRACKET(Ln);

			if ((m_mu*m_epsf > 0) && (m_ss.m_data[i].m_pme))
			{
				// update the metrics
				FESurfaceElement& mel = *m_ss.m_data[i].m_pme;

				double r = m_ss.m_data[i].m_rs[0], s = m_ss.m_data[i].m_rs[1];
				double rp = m_ss.m_data[i].m_rsp[0], sp = m_ss.m_data[i].m_rsp[1];

				double Ln = m_ss.m_data[i].m_Lm;

				mat2d Mk = m_ss.m_data[i].m_M;
				mat2d Mki = Mk.inverse();
				
			
				// update traction multipliers
				// a. trial state
				Lt[0] = m_ss.m_data[i].m_Lt[0] + m_epsf*(Mk[0][0]*(r - rp) + Mk[0][1]*(s - sp));
				Lt[1] = m_ss.m_data[i].m_Lt[1] + m_epsf*(Mk[1][0]*(r - rp) + Mk[1][1]*(s - sp));

				double TMT = Lt[0]*(Mki[0][0]*Lt[0]+Mki[0][1]*Lt[1])+Lt[1]*(Mki[1][0]*Lt[0]+Mki[1][1]*Lt[1]);
				assert(TMT >= 0);

				double phi = sqrt(TMT) - m_mu*Ln;

				// b. return map
				if (phi > 0)
				{
					Lt[0] = m_mu*Ln*Lt[0]/sqrt(TMT);
					Lt[1] = m_mu*Ln*Lt[1]/sqrt(TMT);
				}

				m_ss.m_data[i].m_M = m_ss.Metric0(mel, r, s);

				m_ss.m_data[i].m_Lt[0] = Lt[0];
				m_ss.m_data[i].m_Lt[1] = Lt[1];
			}
		}	

		for (int i=0; i<m_ms.Nodes(); ++i)
		{
			eps = m_ms.m_data[i].m_eps*scale;

			// update Lagrange multipliers
			Ln = m_ms.m_data[i].m_Lm + eps*m_ms.m_data[i].m_gap;
			m_ms.m_data[i].m_Lm = MBRACKET(Ln);

			if ((m_mu*m_epsf > 0) && (m_ms.m_data[i].m_pme))
			{
				// update the metrics
				FESurfaceElement& mel = *m_ms.m_data[i].m_pme;

				double r = m_ms.m_data[i].m_rs[0], s = m_ms.m_data[i].m_rs[1];
				double rp = m_ms.m_data[i].m_rsp[0], sp = m_ms.m_data[i].m_rsp[1];

				double Ln = m_ms.m_data[i].m_Lm;

				mat2d Mk = m_ms.m_data[i].m_M;
				mat2d Mki = Mk.inverse();
				
				// update traction multipliers
				// a. trial state
				double Lt[2];
				Lt[0] = m_ms.m_data[i].m_Lt[0] + m_epsf*(Mk[0][0]*(r - rp) + Mk[0][1]*(s - sp));
				Lt[1] = m_ms.m_data[i].m_Lt[1] + m_epsf*(Mk[1][0]*(r - rp) + Mk[1][1]*(s - sp));

				double TMT = Lt[0]*(Mki[0][0]*Lt[0]+Mki[0][1]*Lt[1])+Lt[1]*(Mki[1][0]*Lt[0]+Mki[1][1]*Lt[1]);
				assert(TMT >= 0);

				double phi = sqrt(TMT) - m_mu*Ln;

				// b. return map
				if (phi > 0)
				{
					Lt[0] = m_mu*Ln*Lt[0]/sqrt(TMT);
					Lt[1] = m_mu*Ln*Lt[1]/sqrt(TMT);
				}

				m_ms.m_data[i].m_M = m_ms.Metric0(mel, r, s);

				m_ms.m_data[i].m_Lt[0] = Lt[0];
				m_ms.m_data[i].m_Lt[1] = Lt[1];
			}
		}
	}

	if (bconv)
	{
		for (int i=0; i<m_ss.Nodes(); ++i) m_ss.m_data[i].m_rsp = m_ss.m_data[i].m_rs;
		for (int i=0; i<m_ms.Nodes(); ++i) m_ms.m_data[i].m_rsp = m_ms.m_data[i].m_rs;
	}

	// store the last gap norm
	m_normg0 = normg1;

	return bconv;
}

//-----------------------------------------------------------------------------
//! This function transforms friction data between two segments

void FESlidingInterface::MapFrictionData(int inode, FESlidingSurface& ss, FESlidingSurface& ms, FESurfaceElement &en, FESurfaceElement &eo, vec3d &q)
{
	// first we find the projection of the old point on the new segment
	double r = ss.m_data[inode].m_rs[0];
	double s = ss.m_data[inode].m_rs[1];
	double rp = ss.m_data[inode].m_rsp[0], ro = rp;
	double sp = ss.m_data[inode].m_rsp[1], so = sp;
	vec3d xn = ms.Local2Global(eo, rp, sp);
	vec3d qn;
	qn = ms.ProjectToSurface(en, xn, rp, sp);
	ss.m_data[inode].m_rsp[0] = rp;
	ss.m_data[inode].m_rsp[1] = sp;

	// next, we transform the frictional traction
	// since these tractions are decomposed in the local 
	// element coordinate system, we have to do a coordinate transformation
	// note that this transformation needs to be done in curvilinear
	// coordinates since the base vectors may not be orthonormal. Also
	// note that we are doing this in the reference configuration
	vec3d to[2], tn[2];
	ms.ContraBaseVectors0(eo, ro, so, to);
	ms.CoBaseVectors0(en, r, s, tn);

	double Lt[2];
	Lt[0] = ss.m_data[inode].m_Lt[0];
	Lt[1] = ss.m_data[inode].m_Lt[1];
	
	vec3d t;
	t = to[0]*Lt[0] + to[1]*Lt[1];

	Lt[0] = t*tn[0];
	Lt[1] = t*tn[1];

	ss.m_data[inode].m_Lt[0] = Lt[0];
	ss.m_data[inode].m_Lt[1] = Lt[1];
}

//-----------------------------------------------------------------------------
void FESlidingInterface::UpdateContactPressures()
{
	int npass = (m_btwo_pass?2:1);
	for (int np=0; np<npass; ++np)
	{
		FESlidingSurface& ss = (np == 0? m_ss : m_ms);
		FESlidingSurface& ms = (np == 0? m_ms : m_ss);
		
		// loop over all nodes of the primary surface
		for (int n=0; n<ss.Nodes(); ++n)
		{
			// get the normal tractions at the integration points
			double gap = ss.m_data[n].m_gap;
			double eps = m_eps*ss.m_data[n].m_eps;
			ss.m_data[n].m_Ln = MBRACKET(ss.m_data[n].m_Lm + eps*gap);
			FESurfaceElement* pme = ss.m_data[n].m_pme;
			if (m_btwo_pass && pme)
			{
				int me = pme->Nodes();
				if (me < 6)
				{
					double ti[6];
					for (int j=0; j<me; ++j) {
						int k = pme->m_lnode[j];
						gap = ms.m_data[k].m_gap;
						eps = m_eps*ms.m_data[k].m_eps;
						ti[j] = MBRACKET(ms.m_data[k].m_Lm + m_eps*ms.m_data[k].m_eps*ms.m_data[k].m_gap);
					}
					// project the data to the nodes
					double tn[6];
					pme->FEElement::project_to_nodes(ti, tn);
					// now evaluate the traction at the intersection point
					double Ln = pme->eval(tn, ss.m_data[n].m_rs[0], ss.m_data[n].m_rs[1]);
					ss.m_data[n].m_Ln += MBRACKET(Ln);
				}
			}
		}
	}
}

//-----------------------------------------------------------------------------
void FESlidingInterface::Serialize(DumpStream& ar)
{
	// store contact data
	FEContactInterface::Serialize(ar);

	// store contact surface data
	m_ms.Serialize(ar);
	m_ss.Serialize(ar);

	// restore element pointers
	SerializePointers(m_ss, m_ms, ar);
	SerializePointers(m_ms, m_ss, ar);

	ar & m_bfirst & m_normg0;
}

void FESlidingInterface::SerializePointers(FESlidingSurface& ss, FESlidingSurface& ms, DumpStream& ar)
{
	if (ar.IsSaving())
	{
		for (int i = 0; i < ss.m_data.size(); i++)
		{
			FESurfaceElement* pe = ss.m_data[i].m_pme;
			int eid = (pe ? pe->m_lid : -1);
			ar << eid;
		}
	}
	else
	{
		for (int i = 0; i < ss.m_data.size(); i++)
		{
			int eid = -1;
			ar >> eid;
			if (eid >= 0)
			{
				FESurfaceElement* pe = &ms.Element(eid);
				ss.m_data[i].m_pme = pe;
			}
			else ss.m_data[i].m_pme = nullptr;
		}
	}
}
