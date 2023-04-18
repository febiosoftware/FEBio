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
#include "FE2OMicroConstraint.h"
#include <FECore/log.h>
#include <FEBioMech/FEBioMech.h>
#include <FECore/FELinearSystem.h>
#include <FECore/FEMesh.h>

//-----------------------------------------------------------------------------
//! constructor
FEMicroFlucSurface::FEMicroFlucSurface(FEModel* fem) : FESurface(fem)
{
	m_Lm.x = 0.; m_Lm.y = 0.; m_Lm.z = 0.;
	m_pv.x = 0.; m_pv.y = 0.; m_pv.z = 0.;
	m_c.x = 0.;  m_c.y = 0.;  m_c.z = 0.;

	m_Fm.unit(); m_Gm.zero();
}

//-----------------------------------------------------------------------------
bool FEMicroFlucSurface::Init()
{
	// calculate the intial microfluctations across the surface
	m_c = SurfMicrofluc();

	return true;
}

//-----------------------------------------------------------------------------
void FEMicroFlucSurface::CopyFrom(FEMicroFlucSurface& s)
{
	m_Node = s.m_Node;

	// create elements
	int NE = s.Elements();
	Create(NE);
	for (int i=0; i<NE; ++i) Element(i) = s.Element(i);

	// copy surface data
	m_Lm = s.m_Lm;
	m_pv = s.m_pv;
	m_c  = s.m_c;
	m_Fm = s.m_Fm;
	m_Gm = s.m_Gm;
}

//-----------------------------------------------------------------------------
//! Calculate the initial volume
vec3d FEMicroFlucSurface::SurfMicrofluc()
{
	// Integration of microfluctation field across surface
	vec3d c;
	
	// get the mesh
	FEMesh& mesh = *GetMesh();

	// loop over all elements
	double vol = 0.0;
	int NE = Elements();
	vec3d x[FEElement::MAX_NODES];
	vec3d x0[FEElement::MAX_NODES];

	for (int i=0; i<NE; ++i)
	{
		// get the next element
		FESurfaceElement& el = Element(i);

		// get the nodal coordinates
		int neln = el.Nodes();
		for (int j=0; j<neln; ++j){
			x[j] = mesh.Node(el.m_node[j]).m_rt;
			x0[j] = mesh.Node(el.m_node[j]).m_r0;
		}

		// loop over integration points
		double* w = el.GaussWeights();
		int nint = el.GaussPoints();
		for (int n=0; n<nint; ++n)
		{
			vec3d r = el.eval(x, n);
			vec3d r0 = el.eval(x0, n);

			vec3d u = r - r0;
			mat3d I; I.unit();

			c += (u - (m_Fm - I)*r0 - m_Gm.contractdyad1(r0)*0.5)*w[n];
		}
	}
	
	return c;
}

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FE2OMicroConstraint, FESurfaceConstraint);
	ADD_PARAMETER(m_blaugon, "laugon" ); 
	ADD_PARAMETER(m_atol   , "augtol" );
	ADD_PARAMETER(m_eps    , "penalty");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! constructor. Set default parameter values
FE2OMicroConstraint::FE2OMicroConstraint(FEModel* pfem) : FESurfaceConstraint(pfem), m_s(pfem), m_dofU(pfem)
{
	m_eps = 0.0;
	m_atol = 0.0;
	m_blaugon = false;
	m_binit = false;	// will be set to true during activation

	// TODO: Can this be done in Init, since there is no error checking
	if (pfem)
	{
		m_dofU.AddVariable(FEBioMech::GetVariableName(FEBioMech::DISPLACEMENT));
	}
}

//-----------------------------------------------------------------------------
void FE2OMicroConstraint::CopyFrom(FENLConstraint* plc)
{
	// cast to a periodic boundary
	FE2OMicroConstraint& mc = dynamic_cast<FE2OMicroConstraint&>(*plc);

	// copy parameters
	GetParameterList() = mc.GetParameterList();

	// copy nodes
	m_s.CopyFrom(mc.m_s);
}

//-----------------------------------------------------------------------------
//! Returns the surface
FESurface* FE2OMicroConstraint::GetSurface()
{
	return &m_s;
}

//-----------------------------------------------------------------------------
//! Initializes data structures. 
void FE2OMicroConstraint::Activate()
{
	// don't forget to call base class
	FENLConstraint::Activate();

	// initialize the surface
	if (m_binit == false) m_s.Init();

	// set flag that initial volume is calculated
	m_binit = true;
}

//-----------------------------------------------------------------------------
void FE2OMicroConstraint::BuildMatrixProfile(FEGlobalMatrix& M)
{
	// We don't do anything here since the connectivity of a surface
	// is implied by the domain to which it is attached.
}

//-----------------------------------------------------------------------------
void FE2OMicroConstraint::UnpackLM(FEElement& el, vector<int>& lm)
{
	FEMesh& mesh = GetMesh();
	int N = el.Nodes();
	lm.resize(N*3);
	for (int i=0; i<N; ++i)
	{
		int n = el.m_node[i];
		FENode& node = mesh.Node(n);
		vector<int>& id = node.m_ID;

		lm[3*i  ] = id[m_dofU[0]];
		lm[3*i+1] = id[m_dofU[1]];
		lm[3*i+2] = id[m_dofU[2]];
	}
}

//-----------------------------------------------------------------------------
void FE2OMicroConstraint::LoadVector(FEGlobalVector& R, const FETimeInfo& tp)
{
	FEMesh& mesh = *m_s.GetMesh();

	vector<double> fe;
	vector<int> lm;

	// get the lagrange 
	vec3d Lm = m_s.m_Lm;

	// loop over all elements
	int NE = m_s.Elements();
	vec3d x[FEElement::MAX_NODES];
	for (int i=0; i<NE; ++i)
	{
		// get the next element
		FESurfaceElement& el = m_s.Element(i);

		// get the nodal coordinates
		int neln = el.Nodes();
		for (int j=0; j<neln; ++j) x[j] = mesh.Node(el.m_node[j]).m_rt;

		// allocate element residual vector
		int ndof = 3*neln;
		fe.resize(ndof);
		zero(fe);

		// loop over all integration points
		double* w = el.GaussWeights();
		int nint = el.GaussPoints();
		for (int n=0; n<nint; ++n)
		{
			// calculate the tangent vectors
			double* Gr = el.Gr(n);
			double* Gs = el.Gs(n);
			vec3d dxr(0,0,0), dxs(0,0,0);
			for (int j=0; j<neln; ++j) 
			{
				dxr += x[j]*Gr[j];
				dxs += x[j]*Gs[j];
			}

			// evaluate the "normal" vector
			vec3d v = (dxr ^ dxs);
			vec3d f = Lm*w[n]*v.norm();

			// evaluate the element forces
			double* H = el.H(n);
			for (int j=0; j<neln; ++j)
			{
				fe[3*j  ] += H[j]*f.x;
				fe[3*j+1] += H[j]*f.y;
				fe[3*j+2] += H[j]*f.z;
			}
		}

		// get the element's LM vector
		lm.resize(3*neln);
		for (int j=0; j<neln; ++j)
		{
			vector<int>& id = mesh.Node(el.m_node[j]).m_ID;
			lm[3*j  ] = id[m_dofU[0]];
			lm[3*j+1] = id[m_dofU[1]];
			lm[3*j+2] = id[m_dofU[2]];
		}

		// add element force vector to global force vector
		R.Assemble(el.m_node, lm, fe);
	}
}

//-----------------------------------------------------------------------------
void FE2OMicroConstraint::StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp)
{
	FEMesh& mesh = *m_s.GetMesh();

	// element stiffness matrix
	FEElementMatrix ke;
	vector<int> lm;
	vector<double> fe;

	// loop over all elements
	int NE = m_s.Elements();
	vec3d x[FEElement::MAX_NODES];
	for (int l=0; l<NE; ++l)
	{
		// get the next element
		FESurfaceElement& el = m_s.Element(l);

		// get the nodal coordinates
		int neln = el.Nodes();
		for (int j=0; j<neln; ++j) x[j] = mesh.Node(el.m_node[j]).m_rt;

		// allocate the stiffness matrix
		int ndof = 3*neln;
		ke.resize(ndof, ndof);
		ke.zero();
		fe.resize(ndof);
		zero(fe);

		// repeat over integration points
		double* w = el.GaussWeights();
		int nint = el.GaussPoints();
		for (int n=0; n<nint; ++n)
		{
			// calculate tangent vectors
			double* N = el.H(n);
			double* Gr = el.Gr(n);
			double* Gs = el.Gs(n);
			vec3d dxr(0,0,0), dxs(0,0,0);
			for (int j=0; j<neln; ++j) 
			{
				dxr += x[j]*Gr[j];
				dxs += x[j]*Gs[j];
			}

			// calculate pressure contribution
			vec3d v = (dxr ^ dxs);

			double vi; 
			double vj;
			for (int i=0; i<neln; ++i)
				for (int j=0; j<neln; ++j)
				{
					vi = N[i]*v.norm();
					vj = N[j]*v.norm();
					ke[3*i  ][3*j  ] += m_eps*vi*vj;
					ke[3*i+1][3*j+1] += m_eps*vi*vj;
					ke[3*i+2][3*j+2] += m_eps*vi*vj;
				}

		
			// calculate displacement contribution
			vec3d qab;
			for (int i=0; i<neln; ++i)
				for (int j=0; j<neln; ++j)
				{
					qab = (-dxs*Gr[j] + dxr*Gs[j])*(N[i]/(2*(dxr ^ dxs).norm()))*w[n]; 

					ke[3*i  ][3*j  ] +=      0;
					ke[3*i  ][3*j+1] +=  qab.z;
					ke[3*i  ][3*j+2] += -qab.y;

					ke[3*i+1][3*j  ] += -qab.z;
					ke[3*i+1][3*j+1] +=      0;
					ke[3*i+1][3*j+2] +=  qab.x;

					ke[3*i+2][3*j  ] +=  qab.y;
					ke[3*i+2][3*j+1] += -qab.x;
					ke[3*i+2][3*j+2] +=      0;
				}
		}


		// get the element's LM vector
		lm.resize(3*neln);
		for (int j=0; j<neln; ++j)
		{
			vector<int>& id = mesh.Node(el.m_node[j]).m_ID;
			lm[3*j  ] = id[m_dofU[0]];
			lm[3*j+1] = id[m_dofU[1]];
			lm[3*j+2] = id[m_dofU[2]];
		}

		// assemble element matrix in global stiffness matrix
		ke.SetNodes(el.m_node);
		ke.SetIndices(lm);
		LS.Assemble(ke);
	}
}

//-----------------------------------------------------------------------------
bool FE2OMicroConstraint::Augment(int naug, const FETimeInfo& tp)
{
	// make sure we are augmenting
	if ((m_blaugon == false) || (m_atol <= 0.0)) return true;

	feLog("\n2O periodic surface microfluctation constraint:\n");

	vec3d Dm = m_s.m_c*m_eps;
	vec3d Lm = m_s.m_pv;
	
	double Dnorm = Dm.norm();
	double Lnorm = Lm.norm();

	double err = Dnorm/Lnorm;

	if (Lnorm == 0)
		err = 0;

	feLog("\tpressure vect norm: %lg\n", Lm.norm());
	feLog("\tnorm : %lg (%lg)\n", err, m_atol);
	feLog("\ttotal microfluc norm: %lg\n", m_s.m_c.norm());

	// check convergence
	if (err < m_atol) return true;

	// update Lagrange multiplier (and pressure variable)
	m_s.m_Lm = Lm;
	m_s.m_pv = Lm + Dm;

	return false;
}

//-----------------------------------------------------------------------------
void FE2OMicroConstraint::Serialize(DumpStream& ar)
{
	ar & m_s.m_Lm;
	ar & m_s.m_pv;
	ar & m_s.m_c;
	ar & m_s.m_Fm;
	ar & m_s.m_Gm;
}

//-----------------------------------------------------------------------------
void FE2OMicroConstraint::Reset()
{
}

//-----------------------------------------------------------------------------
// This function is called when the FE model's state needs to be updated.
void FE2OMicroConstraint::Update(const FETimeInfo& tp)
{
	// calculate the current volume
	m_s.m_c = m_s.SurfMicrofluc();
	
	// update pressure variable
	m_s.m_pv = m_s.m_Lm - m_s.m_c*m_eps;
}
