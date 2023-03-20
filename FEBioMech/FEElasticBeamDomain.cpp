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

FEElasticBeamDomain::FEElasticBeamDomain(FEModel* fem) : FEBeamDomain(fem), FEElasticDomain(fem), m_dofs(fem)
{
	
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
	// TODO:
	return false;
}

//! Get the list of dofs on this domain
const FEDofList& FEElasticBeamDomain::GetDOFList() const
{
	return m_dofs;
}

void FEElasticBeamDomain::SetMaterial(FEMaterial* pm)
{
	m_mat = dynamic_cast<FEElasticBeamMaterial*>(pm);
	assert(m_mat);
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
	double L0 = 1; // TODO: reference length of beam

	// loop over integration points
	int nint = el.GaussPoints();
	int ne = el.Nodes();
	double w = L0; // weight is reference length of beam

	for (int n = 0; n < nint; ++n)
	{
		FEElasticBeamMaterialPoint& mp = dynamic_cast<FEElasticBeamMaterialPoint&>(*el.GetMaterialPoint(n));

		// get stress from beam element
		vec3d t = mp.m_t;	// stress
		vec3d m = mp.m_m;	// moment 

		vec3d G0 = mp.m_G0; // = dphi_0/dS

		mat3da S(G0); // skew-symmetric matrix from Grad (TODO: check sign!)

		// integration point 
		double s = 0.0;

		// shape function (at integration point)
		double H[2] = { 0.5 * (1.0 - s), 0.5 * (1.0 + s) };

		// shape function derivative (at integration point)
		double G[2] = { -0.5, 0.5 }; // TODO: multiply by L0?

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
			for (int j = 0; j < 6; ++j) fe[i * 6 + j] += P[j] * w;
		}
	}
}

//! Calculate global stiffness matrix (only contribution from internal force derivative)
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
	int ne = el.Nodes();

	double L0 = 1; // TODO: reference beam length

	// only a single integration point
	double w = L0;

	// integration point 
	double s = 0.0;

	int nint = el.GaussPoints();
	for (int n = 0; n < nint; ++n)
	{
		FEElasticBeamMaterialPoint& mp = dynamic_cast<FEElasticBeamMaterialPoint&>(*el.GetMaterialPoint(n));

		// get stress from beam element
		vec3d t = mp.m_t;	// stress traction
		vec3d m = mp.m_m;	// moment 

		mat3da St(t);
		mat3da Sm(m);

		// shape function (at integration point)
		double H[2] = { 0.5 * (1.0 - s), 0.5 * (1.0 + s) };

		// shape function derivative (at integration point)
		double G[2] = { -0.5, 0.5 }; // TODO: multiply by L0?

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

				// material stiffness
				matrix Se(6, 6);
				Se = ksi_a * c * ksi_bT;

				// geometrical stiffness
				mat3d P = (t & G0) - I * (t * G0);
				matrix Te(6, 6); Te.zero();
				Te.add(0, 3, (-St * (G[a] * H[b])));
				Te.add(3, 0, (St * (H[a] * G[b])));
				Te.add(3, 3, P * (H[a] * H[b]));

				// assemble into ke
				for (int i = 0; i < 6; ++i)
					for (int j = 0; j < 6; ++j)
					{
						ke[6 * a + i][6 * b + j] += (Se(i, j) + Te(i, j)) * w;
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
	// loop over all integration points
	int nint = el.GaussPoints();
	for (int n = 0; n < nint; ++n)
	{
		// get the material point
		FEElasticBeamMaterialPoint& mp = dynamic_cast<FEElasticBeamMaterialPoint&>(*el.GetMaterialPoint(n));

		// update G0 = dphi0/dS
		mp.m_G0;

		// evaluate the stress
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
