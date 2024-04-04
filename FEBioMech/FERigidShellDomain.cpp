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
#include "FERigidShellDomain.h"
#include "FEBodyForce.h"
#include "FERigidMaterial.h"
#include <FECore/FELinearSystem.h>
#include <FECore/FEModel.h>
#include "FEBioMech.h"

//-----------------------------------------------------------------------------
FERigidShellDomainOld::FERigidShellDomainOld(FEModel* pfem) : FEElasticShellDomainOld(pfem) {}

//-----------------------------------------------------------------------------
// NOTE: Although this function doesn't do anything, we need it since 
//       for rigid shell domains we don't want to call the FEElasticShellDomain::Initialize member.
bool FERigidShellDomainOld::Init()
{
	// just call the base class
	return FEShellDomainOld::Init();
}

//-----------------------------------------------------------------------------
// We need to override it since the base class version will not work for rigid domains.
void FERigidShellDomainOld::Reset()
{
	// nothing here
}

//-----------------------------------------------------------------------------
//! Calculate stiffness contributions for rigid shells.
//! Since rigid elements don't generate stress, we don't need to do
//! anything here.
void FERigidShellDomainOld::StiffnessMatrix(FELinearSystem& LS)
{
	// Caught you looking!
}


//-----------------------------------------------------------------------------
//! calculate residual forces for rigid shells
//!
void FERigidShellDomainOld::InternalForces(FEGlobalVector& R)
{
	// Nothing to do.
}

//-----------------------------------------------------------------------------
//! update stresses for rigid shells.
//!
void FERigidShellDomainOld::Update(const FETimeInfo& tp)
{
	// Nothing to see here. Please move on.
}

//-----------------------------------------------------------------------------
void FERigidShellDomainOld::MassMatrix(FELinearSystem& LS, double scale)
{
	// Only crickets here ... 
}

//-----------------------------------------------------------------------------
void FERigidShellDomainOld::InertialForces(FEGlobalVector& R, std::vector<double>& F)
{
	// chirp, chirp ...
}

//=======================================================================
BEGIN_FECORE_CLASS(FERigidShellDomain, FEShellDomain)
	ADD_PARAMETER(m_h0, "shell_thickness");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FERigidShellDomain::FERigidShellDomain(FEModel* fem) : FEShellDomain(fem), FEElasticDomain(fem), m_dof(fem)
{
	m_h0 = 0.0;
	m_pMat = nullptr;

	// TODO: Can this be done in Init, since there is no error checking
	if (fem)
	{
		m_dof.AddVariable(FEBioMech::GetVariableName(FEBioMech::DISPLACEMENT));
	}
}

bool FERigidShellDomain::Create(int nelems, FE_Element_Spec espec)
{
	m_Elem.resize(nelems);
	for (int i = 0; i < nelems; ++i)
	{
		FEShellElement& el = m_Elem[i];
		el.SetLocalID(i);
		el.SetMeshPartition(this);
	}

	if (espec.etype != FE_ELEM_INVALID_TYPE)
		for (int i = 0; i < nelems; ++i) m_Elem[i].SetType(espec.etype);

	return true;
}

//-----------------------------------------------------------------------------
void FERigidShellDomain::Activate()
{
	// don't need to do anything here
}

//-----------------------------------------------------------------------------
void FERigidShellDomain::AssignDefaultShellThickness()
{
	double h0 = m_h0;
	if (h0 <= 0.0) return;

	for (int j = 0; j < Elements(); ++j)
	{
		FEShellElement& el = Element(j);
		int ne = el.Nodes();
		for (int n = 0; n < ne; ++n) el.m_ht[n] = el.m_h0[n] = h0;
	}
}

void FERigidShellDomain::SetMaterial(FEMaterial* pm)
{
	m_pMat = dynamic_cast<FERigidMaterial*>(pm); assert(m_pMat);
	FEShellDomain::SetMaterial(pm);
}

FEMaterial* FERigidShellDomain::GetMaterial()
{ 
	return m_pMat; 
}

void FERigidShellDomain::Update(const FETimeInfo& tp)
{
	// nothing to do
}

void FERigidShellDomain::Reset()
{
	// nothing to do
}

void FERigidShellDomain::InternalForces(FEGlobalVector& R)
{
	// nothing to do
}

void FERigidShellDomain::BodyForce(FEGlobalVector& R, FEBodyForce& bf)
{
	int NS = (int)m_Elem.size();
#pragma omp parallel for
	for (int i = 0; i < NS; ++i)
	{
		// element force vector
		vector<double> fe;
		vector<int> lm;

		// get the element
		FEShellElement& el = m_Elem[i];

		// create the element force vector and initialize to zero
		int ndof = 3 * el.Nodes();
		fe.assign(ndof, 0);

		// apply body forces to shells
		ElementBodyForce(bf, el, fe);

		// get the element's LM vector
		UnpackLM(el, lm);

		// assemble the residual
		// NOTE: Setting bdom parameter to false to avoid crash, but
		//       need to look further what effect this really has. 
		R.Assemble(el.m_node, lm, fe, false);
	}
}

//-----------------------------------------------------------------------------
//! Calculates element body forces for shells

void FERigidShellDomain::ElementBodyForce(FEBodyForce& BF, FEShellElement& el, vector<double>& fe)
{
	// integration weights
	double* gw = el.GaussWeights();

	// loop over integration points
	int nint = el.GaussPoints();
	int neln = el.Nodes();

	for (int n = 0; n < nint; ++n)
	{
		FEMaterialPoint& mp = *el.GetMaterialPoint(n);
		double dens = m_pMat->Density(mp);

		// calculate the jacobian
		double Jw = detJ0(el, n) * gw[n];

		double* M = el.H(n);

		// get the force
		vec3d f = BF.force(mp);

		for (int i = 0; i < neln; ++i)
		{
			vec3d fi = f * dens * M[i] * Jw;

			fe[3 * i    ] -= fi.x;
			fe[3 * i + 1] -= fi.y;
			fe[3 * i + 2] -= fi.z;
		}
	}
}

void FERigidShellDomain::InertialForces(FEGlobalVector& R, std::vector<double>& F)
{
	// nothing to do
}

void FERigidShellDomain::StiffnessMatrix(FELinearSystem& LS)
{
	// nothing to do
}

void FERigidShellDomain::BodyForceStiffness(FELinearSystem& LS, FEBodyForce& bf)
{
	// repeat over all shell elements
	int NE = (int)m_Elem.size();
#pragma omp parallel for shared (NE)
	for (int iel = 0; iel < NE; ++iel)
	{
		FEShellElement& el = m_Elem[iel];

		// create the element's stiffness matrix
		FEElementMatrix ke(el);
		int ndof = 3 * el.Nodes();
		ke.resize(ndof, ndof);
		ke.zero();

		// calculate inertial stiffness
		ElementBodyForceStiffness(bf, el, ke);

		// get the element's LM vector
		vector<int> lm;
		UnpackLM(el, lm);
		ke.SetIndices(lm);

		// assemble element matrix in global stiffness matrix
		LS.Assemble(ke);
	}
}

//-----------------------------------------------------------------------------
//! This function calculates the stiffness due to body forces
void FERigidShellDomain::ElementBodyForceStiffness(FEBodyForce& BF, FEShellElement &el, matrix &ke)
{
	int neln = el.Nodes();

	double* gw = el.GaussWeights();

	double alphaf = GetFEModel()->GetTime().alphaf;

	// loop over integration points
	int nint = el.GaussPoints();
	for (int n=0; n<nint; ++n)
	{
		FEMaterialPoint& mp = *el.GetMaterialPoint(n);
		double Jw = detJ0(el, n) * gw[n] * alphaf;

		// get the stiffness
		mat3d K = BF.stiffness(mp)*(m_pMat->Density(mp)*Jw);

		double* M = el.H(n);

		for (int i=0, i3=0; i<neln; ++i, i3 += 3)
		{
			for (int j=0, j3 = 0; j<neln; ++j, j3 += 3)
			{
				mat3d Kij = K*(M[i]*M[j]);

				ke[i3  ][j3  ] += Kij(0,0); ke[i3  ][j3+1] += Kij(0,1); ke[i3  ][j3+2] += Kij(0,2);
				ke[i3+1][j3  ] += Kij(1,0); ke[i3+1][j3+1] += Kij(1,1); ke[i3+1][j3+2] += Kij(1,2);
				ke[i3+2][j3  ] += Kij(2,0); ke[i3+2][j3+1] += Kij(2,1); ke[i3+2][j3+2] += Kij(2,2);
			}
		}
	}
}

void FERigidShellDomain::MassMatrix(FELinearSystem& LS, double scale)
{
	// nothing to do here
}

double FERigidShellDomain::detJ0(FEShellElement& el, int n)
{
	vector<vec3d> X(FEElement::MAX_NODES);
	for (int i = 0; i < el.Nodes(); ++i) X[i] = Node(el.m_lnode[i]).m_r0;

	double* Gr = el.Hr(n);
	double* Gs = el.Hs(n);

	vec3d dr = vec3d(0, 0, 0);
	vec3d ds = vec3d(0, 0, 0);
	for (int i = 0; i < el.Nodes(); ++i)
	{
		dr += X[i] * Gr[i];
		ds += X[i] * Gs[i];
	}

	double J0 = (dr ^ ds).norm();

	// multiply with shell thickness
	double h = 0.0;
	double* H = el.H(n);
	for (int i = 0; i < el.Nodes(); ++i) h += H[i] * el.m_h0[i];

	return J0*h;
}

mat3d FERigidShellDomain::CalculateMOI()
{
	mat3d moi; moi.zero();
	FEMesh* pm = GetMesh();

	vector<vec3d> r0(FEElement::MAX_NODES);
	mat3dd I(1); // identity tensor

	// loop over all elements
	for (int iel = 0; iel < Elements(); ++iel)
	{
		FEShellElement& el = Element(iel);

		// initial coordinates
		int neln = el.Nodes();
		for (int i = 0; i < neln; ++i) r0[i] = pm->Node(el.m_node[i]).m_r0;

		// loop over integration points
		double* gw = el.GaussWeights();
		int nint = el.GaussPoints();
		for (int n = 0; n < nint; ++n)
		{
			FEMaterialPoint& mp = *el.GetMaterialPoint(n);

			// calculate jacobian
			double Jw = detJ0(el, n) * gw[n];

			// shape functions at integration point
			double* H = el.H(n);

			double dens = m_pMat->Density(mp);

			// add to moi
			for (int i = 0; i < neln; ++i)
			{
				for (int j = 0; j < neln; ++j)
				{
					mat3d Iij = (r0[i] * r0[j]) * I - (r0[i] & r0[j]);
					moi += Iij * (H[i] * H[j] * Jw * dens);
				}
			}
		}
	}

	return moi;
}

double FERigidShellDomain::CalculateMass()
{
	double mass = 0.0;
	// loop over all elements
	for (int iel = 0; iel < Elements(); ++iel)
	{
		FEShellElement& el = Element(iel);

		// loop over integration points
		int nint = el.GaussPoints();
		double* gw = el.GaussWeights();
		for (int n = 0; n < nint; ++n)
		{
			FEMaterialPoint& mp = *el.GetMaterialPoint(n);

			// calculate jacobian
			double detJ = detJ0(el, n);

			// add to total mass
			mass += m_pMat->Density(mp) * detJ * gw[n];
		}
	}
	return mass;
}

vec3d FERigidShellDomain::CalculateCOM()
{
	vector<vec3d> r0(FEElement::MAX_NODES);
	vec3d rc(0, 0, 0);

	// loop over all elements
	for (int iel = 0; iel < Elements(); ++iel)
	{
		FEShellElement& el = Element(iel);

		// nr of integration points
		int nint = el.GaussPoints();

		// number of nodes
		int neln = el.Nodes();

		// initial coordinates
		for (int i = 0; i < neln; ++i) r0[i] = GetMesh()->Node(el.m_node[i]).m_r0;

		// integration weights
		double* gw = el.GaussWeights();

		// loop over integration points
		for (int n = 0; n < nint; ++n)
		{
			FEMaterialPoint& mp = *el.GetMaterialPoint(n);

			// calculate jacobian
			double detJ = detJ0(el, n);

			// shape functions at integration point
			double* H = el.H(n);

			double dens = m_pMat->Density(mp);

			// add to com
			for (int i = 0; i < el.Nodes(); ++i)
			{
				rc += r0[i] * H[i] * detJ * gw[n] * dens;
			}
		}
	}

	return rc;
}
