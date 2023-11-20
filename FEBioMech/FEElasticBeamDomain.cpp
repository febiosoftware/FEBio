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
#include "FEElasticBeamDomain.h"
#include <FECore/FELinearSystem.h>
#include "FEElasticBeamMaterial.h"
#include "FEBioMech.h"
#include <FECore/FEMesh.h>
#include <FECore/fecore_debug.h>


FEElasticBeamDomain::FEElasticBeamDomain(FEModel* fem) : FEBeamDomain(fem), FEElasticDomain(fem), m_dofs(fem)
{
	m_mat = nullptr;

	if (fem)
	{
		m_dofs.AddVariable(FEBioMech::GetVariableName(FEBioMech::DISPLACEMENT));
		m_dofs.AddVariable(FEBioMech::GetVariableName(FEBioMech::SHELL_ROTATION));
	}
}

// return number of beam elements
int FEElasticBeamDomain::Elements() const
{
	return (int)m_Elem.size();
}

//! return a reference to an element
FEElement& FEElasticBeamDomain::ElementRef(int i) { return m_Elem[i]; }
const FEElement& FEElasticBeamDomain::ElementRef(int i) const { return m_Elem[i]; }

FEBeamElement& FEElasticBeamDomain::Element(int i) { return m_Elem[i]; }

// create function
bool FEElasticBeamDomain::Create(int elements, FE_Element_Spec espec)
{
	m_Elem.resize(elements);
	for (int i = 0; i < elements; ++i)
	{
		FEBeamElement& el = m_Elem[i];
		el.SetLocalID(i);
		el.SetMeshPartition(this);
	}

	// set element type
	int etype = (espec.eshape == ET_LINE3 ? FE_BEAM3G2 : FE_BEAM2G1);
	ForEachElement([=](FEElement& el) { el.SetType(etype); });

	return true;
}

bool FEElasticBeamDomain::Init()
{
	if (FEBeamDomain::Init() == false) return false;

	// initialize elements
	for (int i = 0; i < Elements(); ++i)
	{
		FEBeamElement& el = m_Elem[i];
		vec3d r0[FEElement::MAX_NODES];
		int ne = el.Nodes();
		for (int j = 0; j < ne; ++j) r0[j] = Node(el.m_lnode[j]).m_r0;

		// NOTE: This assumes the beam is initially straight!
		//       This also assumes that nodes 0 and 1 define the boundary nodes. 
		el.m_L0 = (r0[1] - r0[0]).Length();

		vec3d e = r0[1] - r0[0]; e.Normalize();
		el.m_E3 = e;
	}

	return true;
}

//! Get the list of dofs on this domain
const FEDofList& FEElasticBeamDomain::GetDOFList() const
{
	return m_dofs;
}

void FEElasticBeamDomain::SetMaterial(FEMaterial* pm)
{
	FEDomain::SetMaterial(pm);
	m_mat = dynamic_cast<FEElasticBeamMaterial*>(pm);
	assert(m_mat);
}

FEMaterial* FEElasticBeamDomain::GetMaterial()
{
	return m_mat;
}

void FEElasticBeamDomain::PreSolveUpdate(const FETimeInfo& tp)
{
	for (int iel = 0; iel < Elements(); ++iel)
	{
		FEBeamElement& el = Element(iel);
		if (el.isActive())
		{
			int nint = el.GaussPoints();
			for (int n = 0; n < nint; ++n) el.GetMaterialPoint(n)->Update(tp);
		}
	}
}

//! calculate the internal forces
void FEElasticBeamDomain::InternalForces(FEGlobalVector& R)
{
	for (int iel = 0; iel < Elements(); ++iel)
	{
		FEBeamElement& el = Element(iel);

		int ne = el.Nodes();
		int ndof = ne * 6;
		vector<double> fe(ndof, 0.0);

		ElementInternalForces(el, fe);

		vector<int> lm;
		UnpackLM(el, lm);
		R.Assemble(lm, fe);
	}
}

void FEElasticBeamDomain::ElementInternalForces(FEBeamElement& el, std::vector<double>& fe)
{
	// reference length of beam
	double L0 = el.m_L0;

	// loop over integration points
	int nint = el.GaussPoints();
	int ne = el.Nodes();
	double* w = el.GaussWeights();
	for (int n = 0; n < nint; ++n)
	{
		FEElasticBeamMaterialPoint& mp = *(el.GetMaterialPoint(n)->ExtractData<FEElasticBeamMaterialPoint>());

		// shape functions
		double* H = el.H(n);
		double* Hr = el.Hr(n);

		// get stress from beam element
		vec3d t = mp.m_t;	// stress
		vec3d m = mp.m_m;	// moment 

		vec3d G0 = mp.m_G0; // = dphi_0/dS

		mat3da S(G0); // skew-symmetric matrix from Grad (TODO: check sign!)

		// J = dS / dr
		double J = L0 / 2.0;

		// shape function derivative (at integration point)
		double G[FEElement::MAX_NODES];
		for (int i = 0; i < ne; ++i) G[i] = Hr[i] / J;

		double wJ = w[n] * J;

		for (int i = 0; i < ne; ++i)
		{
			// build ksi matrix
			mat3dd I(1.0);
			matrix ksi(6, 6); ksi.zero();
			ksi.add(0, 0, (I * G[i]));
			ksi.add(3, 3, (I * G[i]));
			ksi.add(3, 0, (-S * H[i]));

			vector<double> R{ t.x, t.y, t.z, m.x, m.y, m.z };
			vector<double> P(6);
			P = ksi * R;

			// assemble
			// note the negative sign: this is because we need to subtract
			// the internal forces from the external forces.
			for (int j = 0; j < 6; ++j) fe[i * 6 + j] -= P[j] * wJ;
		}
	}
}

//! Calculate global stiffness matrix
void FEElasticBeamDomain::StiffnessMatrix(FELinearSystem& LS)
{
	for (int iel = 0; iel < Elements(); ++iel)
	{
		FEBeamElement& el = Element(iel);

		int ne = el.Nodes();
		int ndof = ne * 6;
		vector<int> lm;
		UnpackLM(el, lm);

		FEElementMatrix ke(el, lm); ke.zero();
		ElementStiffnessMatrix(el, ke);

		LS.Assemble(ke);
	}
}

void FEElasticBeamDomain::ElementStiffnessMatrix(FEBeamElement& el, FEElementMatrix& ke)
{
	// reference length of beam
	double L0 = el.m_L0;

	// loop over integration points
	int nint = el.GaussPoints();
	int ne = el.Nodes();
	double* w = el.GaussWeights();
	for (int n = 0; n < nint; ++n)
	{
		FEElasticBeamMaterialPoint& mp = *(el.GetMaterialPoint(n)->ExtractData<FEElasticBeamMaterialPoint>());

		// get stress from beam element
		vec3d t = mp.m_t;	// stress traction
		vec3d m = mp.m_m;	// moment 

		mat3da St(t);
		mat3da Sm(m);

		// shape functions
		double* H = el.H(n);
		double* Hr = el.Hr(n);

		// J = dS / dr
		double J = L0 / 2.0;

		// shape function derivative (at integration point)
		double G[FEElement::MAX_NODES];
		for (int i = 0; i < ne; ++i) G[i] = Hr[i] / J;

		double wJ = w[n] * J;

		vec3d G0 = mp.m_G0; // = dphi_0/dS

		mat3da S(G0); // skew-symmetric matrix from Grad (TODO: check sign!)

		for (int a = 0; a < ne; ++a)
			for (int b = 0; b < ne; ++b)
			{
				// build ksi matrix
				mat3dd I(1.0);
				matrix ksi_a(6, 6); ksi_a.zero();
				matrix ksi_bT(6, 6); ksi_bT.zero();
				ksi_a.add(0, 0, (I * G[a]));
				ksi_a.add(3, 3, (I * G[a]));
				ksi_a.add(3, 0, (-S * H[a]));

				// we assemble the transpose for b
				ksi_bT.add(0, 0, (I * G[b]));
				ksi_bT.add(3, 3, (I * G[b]));
				ksi_bT.add(0, 3, (S * H[b]));

				// spatial tangent
				matrix c(6, 6);
				m_mat->Tangent(mp, c);

				// material stiffness
				matrix Se(6, 6);
				Se = ksi_a * c * ksi_bT;

				// geometrical stiffness
				mat3d P = (t & G0) - I * (t * G0);
				matrix Te(6, 6); Te.zero();
				Te.add(0, 3, (-St * (G[a] * H[b])));
				Te.add(3, 0, (St * (H[a] * G[b])));
				Te.add(3, 3, P * (H[a] * H[b]) - Sm*(G[a]*H[b]) );

				// assemble into ke
				for (int i = 0; i < 6; ++i)
					for (int j = 0; j < 6; ++j)
					{
						ke[6 * a + i][6 * b + j] += (Se(i, j) + Te(i, j)) * wJ;
					}
			}
	}
}

void FEElasticBeamDomain::IncrementalUpdate(std::vector<double>& ui, bool finalFlag)
{
	// update the rotations at the material points
	for (int i = 0; i < Elements(); ++i)
	{
		FEBeamElement& el = Element(i);
		int ne = el.Nodes();
		int ni = el.GaussPoints();

		// initial length
		double L0 = el.m_L0;
		double J = L0 / 2.0;

		// get the nodal values of the incremental rotations
		vector<vec3d> dri(ne);
		for (int j = 0; j < ne; ++j)
		{
			std::vector<int>& id = Node(el.m_lnode[j]).m_ID;
			vec3d d(0, 0, 0);
			if (id[m_dofs[3]] >= 0) d.x = ui[id[m_dofs[3]]];
			if (id[m_dofs[4]] >= 0) d.y = ui[id[m_dofs[4]]];
			if (id[m_dofs[5]] >= 0) d.z = ui[id[m_dofs[5]]];

			dri[j] = d;
		}

		// evaluate at integration points
		for (int n = 0; n < ni; ++n)
		{
			// get the material point
			FEElasticBeamMaterialPoint& mp = *(el.GetMaterialPoint(n)->ExtractData<FEElasticBeamMaterialPoint>());

			// evaluate at integration point
			vec3d dr = el.Evaluate(&dri[0], n);

			// use shape function derivatives
			double* Hr = el.Hr(n);
			vec3d drdS(0, 0, 0);
			for (int a = 0; a < ne; ++a) drdS += dri[a] * (Hr[a] / J);

			// calculate exponential map
			mat3d dR; dR.exp(dr);

			// extract quaternion
			quatd dq(dR);

			// update rotations
			mp.m_Rt = (dq*mp.m_Ri) * mp.m_Rp;

			// update spatial curvature
			mat3da Wn(mp.m_wn);
			mat3d W = dR * Wn * dR.transpose(); // this should be a skew-symmetric matrix!

			vec3d w2(-W[1][2], W[0][2], -W[0][1]);

			double g1 = 1, g2 = 1;
			double a = dr.norm();
			if (a != 0)
			{
				g1 = sin(a) / a;
				g2 = sin(0.5 * a) / (0.5 * a);
			}
			vec3d e(dr); e.unit();

			vec3d w1 = drdS*g1 + e*((1.0 - g1)*(e*drdS)) + (dr ^ drdS)*(0.5*g2*g2);

			// update curvature
			mp.m_w = w1 + w2;

			if (finalFlag)
			{
				mp.m_Ri = dq * mp.m_Ri;
				mp.m_wn = mp.m_w;
			}
		}
	}
}

void FEElasticBeamDomain::Update(const FETimeInfo& tp)
{
	int NE = Elements();
	for (int i = 0; i < NE; ++i)
	{
		FEBeamElement& el = Element(i);
		UpdateElement(el);
	}
}

void FEElasticBeamDomain::UpdateElement(FEBeamElement& el)
{
	// get the nodal positions
	vec3d rt[FEElement::MAX_NODES];
	int ne = el.Nodes();
	for (int i = 0; i < ne; ++i) rt[i] = Node(el.m_lnode[i]).m_rt;

	// initial length
	double L0 = el.m_L0;
	double J = L0 / 2;

	// loop over all integration points
	int nint = el.GaussPoints();
	for (int n = 0; n < nint; ++n)
	{
		// get the material point
		FEElasticBeamMaterialPoint& mp = *(el.GetMaterialPoint(n)->ExtractData<FEElasticBeamMaterialPoint>());

		// evaluate G0 = dphi0/dS
		double* Hr = el.Hr(n);
		vec3d G0(0, 0, 0);
		for (int a = 0; a < ne; ++a)
		{
			G0 += rt[a] * (Hr[a] / J);
		}

		// update G0 = dphi0/dS
		mp.m_G0 = G0;
		quatd q = mp.m_Rt;
		quatd qi = q.Conjugate();

		// calculate material strain measures
		mp.m_Gamma = qi * G0 - el.m_E3;
		mp.m_Kappa = qi * mp.m_w; // m_w is updated in IncrementalUpdate(std::vector<double>& ui)

		// evaluate the (spatial) stress
		m_mat->Stress(mp);
	}
}

//! TODO: calculate the interial forces (for dynamic problems)
void FEElasticBeamDomain::InertialForces(FEGlobalVector& R, std::vector<double>& F) { assert(false); }

//! TODO: calculate the mass matrix (for dynamic problems)
void FEElasticBeamDomain::MassMatrix(FELinearSystem& LS, double scale) { assert(false); }

//! TODO: Calculate the body force vector
void FEElasticBeamDomain::BodyForce(FEGlobalVector& R, FEBodyForce& bf) { assert(false); }

//! TODO: Calculate stiffness contribution of body forces
void FEElasticBeamDomain::BodyForceStiffness(FELinearSystem& LS, FEBodyForce& bf) { assert(false); }
