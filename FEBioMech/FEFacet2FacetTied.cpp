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
#include "FEFacet2FacetTied.h"
#include "FECore/FEMesh.h"
#include "FECore/FEClosestPointProjection.h"
#include "FECore/FEGlobalMatrix.h"
#include <FECore/FELinearSystem.h>
#include "FECore/log.h"

//-----------------------------------------------------------------------------
// Define tied interface parameters
BEGIN_FECORE_CLASS(FEFacet2FacetTied, FEContactInterface)
	ADD_PARAMETER(m_atol   , "tolerance"       );
	ADD_PARAMETER(m_eps    , "penalty"         );
	ADD_PARAMETER(m_naugmin, "minaug"          );
	ADD_PARAMETER(m_naugmax, "maxaug"          );
	ADD_PARAMETER(m_stol   , "search_tolerance");
	ADD_PARAMETER(m_gapoff , "gap_offset"      );
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEFacetTiedSurface::Data::Data()
{
	m_vgap  = vec3d(0,0,0);
	m_vgap0 = vec3d(0,0,0);
	m_Lm = vec3d(0,0,0);
	m_rs = vec2d(0,0);
	m_pme = (FESurfaceElement*) 0;
}

void FEFacetTiedSurface::Data::Serialize(DumpStream& ar)
{
	FEContactMaterialPoint::Serialize(ar);
	ar & m_vgap & m_vgap0 & m_Lm & m_rs;
}

//-----------------------------------------------------------------------------
FEFacetTiedSurface::FEFacetTiedSurface(FEModel* pfem) : FEContactSurface(pfem)
{

}

//-----------------------------------------------------------------------------
bool FEFacetTiedSurface::Init()
{
	// initialize surface data first
	if (FEContactSurface::Init() == false) return false;

	return true;
}

//-----------------------------------------------------------------------------
//! create material point data
FEMaterialPoint* FEFacetTiedSurface::CreateMaterialPoint()
{
	return new FEFacetTiedSurface::Data;
}

//-----------------------------------------------------------------------------
void FEFacetTiedSurface::Serialize(DumpStream &ar)
{
	FEContactSurface::Serialize(ar);
	if (ar.IsShallow())
	{
		if (ar.IsSaving())
		{
			for (int n=0; n<Elements(); ++n)
			{
				FESurfaceElement& el = Element(n);
				int nint = el.GaussPoints();
				for (int j=0; j<nint; ++j)
				{
					Data& d = static_cast<Data&>(*el.GetMaterialPoint(j));
					ar << d.m_vgap;
					ar << d.m_vgap0;
					ar << d.m_rs;
					ar << d.m_Lm;
				}
			}
		}
		else
		{
			for (int n = 0; n<Elements(); ++n)
			{
				FESurfaceElement& el = Element(n);
				int nint = el.GaussPoints();
				for (int j = 0; j<nint; ++j)
				{
					Data& d = static_cast<Data&>(*el.GetMaterialPoint(j));
					ar >> d.m_vgap;
					ar >> d.m_vgap0;
					ar >> d.m_rs;
					ar >> d.m_Lm;
				}
			}
		}
	}
	else
	{
		if (ar.IsSaving())
		{
			for (int n = 0; n<Elements(); ++n)
			{
				FESurfaceElement& el = Element(n);
				int nint = el.GaussPoints();
				for (int j = 0; j<nint; ++j)
				{
					Data& d = static_cast<Data&>(*el.GetMaterialPoint(j));
					ar << d.m_vgap;
					ar << d.m_vgap0;
					ar << d.m_rs;
					ar << d.m_Lm;
				}
			}
		}
		else
		{
			for (int n = 0; n<Elements(); ++n)
			{
				FESurfaceElement& el = Element(n);
				int nint = el.GaussPoints();
				for (int j = 0; j<nint; ++j)
				{
					Data& d = static_cast<Data&>(*el.GetMaterialPoint(j));
					ar >> d.m_vgap;
					ar >> d.m_vgap0;
					ar >> d.m_rs;
					ar >> d.m_Lm;
				}
			}
		}
	}
}

//=============================================================================
FEFacet2FacetTied::FEFacet2FacetTied(FEModel* pfem) : FEContactInterface(pfem), m_ss(pfem), m_ms(pfem)
{
	// give this interface an ID
	static int count = 1;
	SetID(count++);

	// set parents
	m_ss.SetContactInterface(this);
	m_ms.SetContactInterface(this);

	// define sibling relationships
	m_ss.SetSibling(&m_ms);
	m_ms.SetSibling(&m_ss);

	// initial parameter values
	m_atol    = 0.01;
	m_eps     = 1.0;
	m_naugmin = 0;
	m_naugmax = 10;
	m_stol    = 0.0001;
	m_gapoff  = false;
}

//-----------------------------------------------------------------------------
//! build the matrix profile for use in the stiffness matrix
void FEFacet2FacetTied::BuildMatrixProfile(FEGlobalMatrix& K)
{
	FEMesh& mesh = GetMesh();

	// get the DOFS
	const int dof_X = GetDOFIndex("x");
	const int dof_Y = GetDOFIndex("y");
	const int dof_Z = GetDOFIndex("z");
	const int dof_RU = GetDOFIndex("Ru");
	const int dof_RV = GetDOFIndex("Rv");
	const int dof_RW = GetDOFIndex("Rw");

	vector<int> lm(6*FEElement::MAX_NODES*2);

	FEFacetTiedSurface& ss = m_ss;
	FEFacetTiedSurface& ms = m_ms;

	for (int j=0; j<ss.Elements(); ++j)
	{
		FESurfaceElement& se = ss.Element(j);
		int nint = se.GaussPoints();
		int* sn = &se.m_node[0];
		for (int k=0; k<nint; ++k)
		{
			FEFacetTiedSurface::Data& d = static_cast<FEFacetTiedSurface::Data&>(*se.GetMaterialPoint(k));
			FESurfaceElement* pe = d.m_pme;
			if (pe != 0)
			{
				FESurfaceElement& me = *pe;
				int* mn = &me.m_node[0];

				lm.assign(lm.size(), -1);

				int nseln = se.Nodes();
				int nmeln = me.Nodes();

				for (int l=0; l<nseln; ++l)
				{
					vector<int>& id = mesh.Node(sn[l]).m_ID;
					lm[6*l  ] = id[dof_X];
					lm[6*l+1] = id[dof_Y];
					lm[6*l+2] = id[dof_Z];
					lm[6*l+3] = id[dof_RU];
					lm[6*l+4] = id[dof_RV];
					lm[6*l+5] = id[dof_RW];
				}

				for (int l=0; l<nmeln; ++l)
				{
					vector<int>& id = mesh.Node(mn[l]).m_ID;
					lm[6*(l+nseln)  ] = id[dof_X];
					lm[6*(l+nseln)+1] = id[dof_Y];
					lm[6*(l+nseln)+2] = id[dof_Z];
					lm[6*(l+nseln)+3] = id[dof_RU];
					lm[6*(l+nseln)+4] = id[dof_RV];
					lm[6*(l+nseln)+5] = id[dof_RW];
				}

				K.build_add(lm);
			}
		}
	}
}

//-----------------------------------------------------------------------------
//! Initialization. This function intializes the surfaces data
bool FEFacet2FacetTied::Init()
{
	// create the surfaces
	if (m_ss.Init() == false) return false;
	if (m_ms.Init() == false) return false;

	return true;
}

//-----------------------------------------------------------------------------
//! Interface activation. Also projects primary surface onto secondary surface
void FEFacet2FacetTied::Activate()
{
	// Don't forget to call base member!
	FEContactInterface::Activate();

	// project primary surface onto secondary surface
	ProjectSurface(m_ss, m_ms);

	if (m_gapoff)
	{
		// loop over all primary elements
		for (int i = 0; i < m_ss.Elements(); ++i)
		{
			// get the primary element
			FESurfaceElement& se = m_ss.Element(i);

			// loop over all its integration points
			int nint = se.GaussPoints();
			for (int j = 0; j < nint; ++j)
			{
				// get integration point data
				FEFacetTiedSurface::Data& pt = static_cast<FEFacetTiedSurface::Data&>(*se.GetMaterialPoint(j));

				pt.m_vgap0 = pt.m_vgap;
				pt.m_vgap = vec3d(0, 0, 0);
				pt.m_gap = 0.0;
			}
		}
	}
}

//-----------------------------------------------------------------------------
void FEFacet2FacetTied::ProjectSurface(FEFacetTiedSurface& ss, FEFacetTiedSurface& ms)
{
	// get the mesh
	FEMesh& mesh = *ss.GetMesh();

	// closest point projection method
	FEClosestPointProjection cpp(ms);
	cpp.HandleSpecialCases(true);
	cpp.SetTolerance(m_stol);
	cpp.Init();

	// loop over all primary elements
	for (int i=0; i<ss.Elements(); ++i)
	{
		// get the primary element
		FESurfaceElement& se = ss.Element(i);

		// get nodal coordinates
		int nn = se.Nodes();
		vec3d re[FEElement::MAX_NODES];
		for (int j=0; j<nn; ++j) re[j] = mesh.Node(se.m_node[j]).m_rt;

		// loop over all its integration points
		int nint = se.GaussPoints();
		for (int j=0; j<nint; ++j)
		{
			// get integration point data
			FEFacetTiedSurface::Data& pt = static_cast<FEFacetTiedSurface::Data&>(*se.GetMaterialPoint(j));

			// calculate the global coordinates of this integration point
			vec3d x = se.eval(re, j);

			// find the secondary element
			vec3d q; vec2d rs;
			FESurfaceElement* pme = cpp.Project(x, q, rs);
			if (pme)
			{
				// store the secondary element
				pt.m_pme = pme;
				pt.m_rs[0] = rs[0];
				pt.m_rs[1] = rs[1];

				// calculate gap
				pt.m_vgap = x - q;

				pt.m_gap = pt.m_vgap.norm();
			}
			else pt.m_pme = 0;
		}
	}
}

//-----------------------------------------------------------------------------
//! Update tied interface data. This function re-evaluates the gaps between
//! the primary node and their projections onto the secondary surface.
//!
void FEFacet2FacetTied::Update()
{
	// get the mesh
	FEMesh& mesh = *m_ss.GetMesh();

	// loop over all primary elements
	const int NE = m_ss.Elements();
	for (int i=0; i<NE; ++i)
	{
		// next element
		FESurfaceElement& se = m_ss.Element(i);
		int nseln = se.Nodes();

		// get the nodal coordinates
		vec3d rs[FEElement::MAX_NODES];
		for (int j=0; j<nseln; ++j) rs[j] = mesh.Node(se.m_node[j]).m_rt;

		// loop over all integration points
		const int nint = se.GaussPoints();
		for (int n=0; n<nint; ++n)
		{
			// get integration point data
			FEFacetTiedSurface::Data& pt = static_cast<FEFacetTiedSurface::Data&>(*se.GetMaterialPoint(n));

			FESurfaceElement* pme = pt.m_pme;
			if (pme)
			{
				FESurfaceElement& me = static_cast<FESurfaceElement&>(*pme);

				// get the current primary nodal position
				vec3d rn = se.eval(rs, n);

				// get the natural coordinates of the primary projection
				double r = pt.m_rs[0];
				double s = pt.m_rs[1];

				// get the secondary nodal coordinates
				int nmeln = me.Nodes();
				vec3d y[FEElement::MAX_NODES];
				for (int l=0; l<nmeln; ++l) y[l] = mesh.Node( me.m_node[l] ).m_rt;

				// calculate the primary node projection
				vec3d q = me.eval(y, r, s);

				// calculate the gap function
				pt.m_vgap = (rn - q) - pt.m_vgap0;

				pt.m_gap = pt.m_vgap.norm();
			}
		}
	}
}


//-----------------------------------------------------------------------------
//! This function calculates the contact forces for a tied interface.
void FEFacet2FacetTied::LoadVector(FEGlobalVector& R, const FETimeInfo& tp)
{
	vector<int> sLM, mLM, LM, en;
	vector<double> fe;

	// shape functions
	double Hm[FEElement::MAX_NODES];

	// get the mesh
	FEMesh& mesh = *m_ss.GetMesh();

	// loop over all elements
	const int NE = m_ss.Elements();
	for (int i=0; i<NE; ++i)
	{
		// get the next element
		FESurfaceElement& se = m_ss.Element(i);
		int nseln = se.Nodes();

		// integration weights
		double* w = se.GaussWeights();

		// get the element's LM vector
		m_ss.UnpackLM(se, sLM);

		// loop over integration points
		const int nint = se.GaussPoints();
		for (int n=0; n<nint; ++n)
		{
			// get integration point data
			FEFacetTiedSurface::Data& pt = static_cast<FEFacetTiedSurface::Data&>(*se.GetMaterialPoint(n));

			// get the secondary element
			FESurfaceElement* pme = pt.m_pme;
			if (pme)
			{
				// get the secondary element
				FESurfaceElement& me = *pme;
				m_ms.UnpackLM(me, mLM);
				int nmeln = me.Nodes();

				// get primary contact force
				vec3d tc = pt.m_Lm + pt.m_vgap*m_eps;

				// calculate jacobian
				// note that we are integrating over the reference surface
				double detJ = m_ss.jac0(se, n);

				// primary shape functions
				double* Hs = se.H(n);

				// secondary shape functions
				double r = pt.m_rs[0];
				double s = pt.m_rs[1];
				me.shape_fnc(Hm, r, s);

				// calculate degrees of freedom
				int ndof = 3*(nseln + nmeln);

				// calculate the force vector
				fe.resize(ndof);
				for (int k=0; k<nseln; ++k)
				{
					fe[3*k  ] = -detJ*w[n]*tc.x*Hs[k];
					fe[3*k+1] = -detJ*w[n]*tc.y*Hs[k];
					fe[3*k+2] = -detJ*w[n]*tc.z*Hs[k];
				}
				for (int k=0; k<nmeln; ++k)
				{
					fe[3*(k+nseln)  ] = detJ*w[n]*tc.x*Hm[k];
					fe[3*(k+nseln)+1] = detJ*w[n]*tc.y*Hm[k];
					fe[3*(k+nseln)+2] = detJ*w[n]*tc.z*Hm[k];
				}

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

				// assemble the global residual
				R.Assemble(en, LM, fe);
			}
		}
	}
}

//-----------------------------------------------------------------------------
//! Calculate the stiffness matrix contribution.
void FEFacet2FacetTied::StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp)
{
	vector<int> sLM, mLM, LM, en;
	FEElementMatrix ke;

	// shape functions
	double Hm[FEElement::MAX_NODES];

	// loop over all primary elements
	const int NE = m_ss.Elements();
	for (int i=0; i<NE; ++i)
	{
		// get the next element
		FESurfaceElement& se = m_ss.Element(i);
		int nseln = se.Nodes();

		// get the element's LM vector
		m_ss.UnpackLM(se, sLM);

		// integration weights
		double* w = se.GaussWeights();

		// loop over all integration points
		const int nint = se.GaussPoints();
		for (int n=0; n<nint; ++n)
		{
			// get intgration point data
			FEFacetTiedSurface::Data& pt = static_cast<FEFacetTiedSurface::Data&>(*se.GetMaterialPoint(n));

			// get the secondary element
			FESurfaceElement* pme = pt.m_pme;
			if (pme)
			{
				// get the secondary element
				FESurfaceElement& me = *pme;
				int nmeln = me.Nodes();
				m_ms.UnpackLM(me, mLM);

				// calculate jacobian
				double detJ = m_ss.jac0(se, n);

				// primary shape functions
				double* Hs = se.H(n);

				// secondary shape functions
				double r = pt.m_rs[0];
				double s = pt.m_rs[1];
				me.shape_fnc(Hm, r, s);

				// calculate degrees of freedom
				int ndof = 3*(nseln + nmeln);

				// create the stiffness matrix
				ke.resize(ndof, ndof);
				ke.zero();
				for (int k=0; k<nseln; ++k)
				{
					for (int l=0; l<nseln; ++l)
					{
						ke[3*k  ][3*l  ] = Hs[k]*Hs[l];
						ke[3*k+1][3*l+1] = Hs[k]*Hs[l];
						ke[3*k+2][3*l+2] = Hs[k]*Hs[l];
					}
				}

				for (int k=0; k<nseln; ++k)
				{
					for (int l=0; l<nmeln; ++l)
					{
						ke[3*k  ][3*(l+nseln)  ] = -Hs[k]*Hm[l];
						ke[3*k+1][3*(l+nseln)+1] = -Hs[k]*Hm[l];
						ke[3*k+2][3*(l+nseln)+2] = -Hs[k]*Hm[l];

						ke[3*(l+nseln)  ][3*k  ] = -Hs[k]*Hm[l];
						ke[3*(l+nseln)+1][3*k+1] = -Hs[k]*Hm[l];
						ke[3*(l+nseln)+2][3*k+2] = -Hs[k]*Hm[l];
					}
				}

				for (int k=0; k<nmeln; ++k)
					for (int l=0; l<nmeln; ++l)
					{
						ke[3*(k+nseln)  ][3*(l+nseln)  ] = Hm[k]*Hm[l];
						ke[3*(k+nseln)+1][3*(l+nseln)+1] = Hm[k]*Hm[l];
						ke[3*(k+nseln)+2][3*(l+nseln)+2] = Hm[k]*Hm[l];
					}

				for (int k=0; k<ndof; ++k)
					for (int l=0; l<ndof; ++l) ke[k][l] *= m_eps*detJ*w[n];

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

				// assemble the global residual
				ke.SetNodes(en);
				ke.SetIndices(LM);
				LS.Assemble(ke);
			}
		}
	}
}

//-----------------------------------------------------------------------------
//! Do an augmentation.
bool FEFacet2FacetTied::Augment(int naug, const FETimeInfo& tp)
{
	// make sure we need to augment
	if (m_laugon != 1) return true;

	// calculate initial norms
	double normL0 = 0;
	for (int i=0; i<m_ss.Elements(); ++i)
	{
		FESurfaceElement& se = m_ss.Element(i);
		for (int j=0; j<se.GaussPoints(); ++j)
		{
			FEFacetTiedSurface::Data& pt = static_cast<FEFacetTiedSurface::Data&>(*se.GetMaterialPoint(j));
			vec3d& lm = pt.m_Lm;
			normL0 += lm*lm;
		}
	}
	normL0 = sqrt(normL0);

	// update Lagrange multipliers and calculate current norms
	double normL1 = 0;
	double normgc = 0;
	int N = 0;
	for (int i=0; i<m_ss.Elements(); ++i)
	{
		FESurfaceElement& se = m_ss.Element(i);
		for (int j=0; j<se.GaussPoints(); ++j)
		{
			FEFacetTiedSurface::Data& pt = static_cast<FEFacetTiedSurface::Data&>(*se.GetMaterialPoint(j));

			vec3d lm = pt.m_Lm + pt.m_vgap*m_eps;

			normL1 += lm*lm;
			if (pt.m_pme != 0)
			{
				double g = pt.m_vgap.norm();
				normgc += g*g;
				++N;
			}
		}
	}	
	if (N == 0) N=1;
	normL1 = sqrt(normL1);
	normgc = sqrt(normgc / N);

	// check convergence of constraints
	feLog(" tied interface # %d\n", GetID());
	feLog("                        CURRENT        REQUIRED\n");
	double pctn = 0;
	if (fabs(normL1) > 1e-10) pctn = fabs((normL1 - normL0)/normL1);
	feLog("    normal force : %15le %15le\n", pctn, m_atol);
	feLog("    gap function : %15le       ***\n", normgc);
		
	// check convergence
	bool bconv = true;
	if (pctn >= m_atol) bconv = false;
	if (naug < m_naugmin ) bconv = false;
	if (naug >= m_naugmax) bconv = true;

	if (bconv == false) 
	{
		for (int i=0; i<m_ss.Elements(); ++i)
		{
			FESurfaceElement& se = m_ss.Element(i);
			for (int j=0; j<se.GaussPoints(); ++j)
			{
				FEFacetTiedSurface::Data& pt = static_cast<FEFacetTiedSurface::Data&>(*se.GetMaterialPoint(j));

				// update Lagrange multipliers
				pt.m_Lm = pt.m_Lm + pt.m_vgap*m_eps;
			}
		}	
	}

	return bconv;
}

//-----------------------------------------------------------------------------
//! Serialize the data to the archive.
void FEFacet2FacetTied::Serialize(DumpStream &ar)
{
	// store contact data
	FEContactInterface::Serialize(ar);

	// store contact surface data
	m_ss.Serialize(ar);
	m_ms.Serialize(ar);
}
