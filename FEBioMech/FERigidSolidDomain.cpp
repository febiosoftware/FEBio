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
#include "FERigidSolidDomain.h"
#include "FERigidMaterial.h"
#include <FECore/FEMesh.h>
#include "FEBodyForce.h"

//-----------------------------------------------------------------------------
FERigidSolidDomain::FERigidSolidDomain(FEModel* pfem) : FEElasticSolidDomain(pfem) {}

//-----------------------------------------------------------------------------
// NOTE: Although this function doesn't do anything we need it because we don't
//       want to call the FEElasticSolidDomain::Initialize function.
bool FERigidSolidDomain::Init()
{
	if (FESolidDomain::Init() == false) return false;

	// set the rigid nodes ID
	// Altough this is already done in FERigidSystem::CreateObjects()
	// we do it here again since the mesh adaptors will call this function
	// and they may have added some nodes
	FERigidMaterial* pm = dynamic_cast<FERigidMaterial*>(m_pMat); assert(pm);
	if (pm == nullptr) return false;
	int rbid = pm->GetRigidBodyID();
	for (int i = 0; i < Nodes(); ++i)
	{
		FENode& node = Node(i);
		node.m_rid = rbid;
	}
	return true;
}

//-----------------------------------------------------------------------------
// We need to override it since the base class version will not work for rigid domains.
void FERigidSolidDomain::Reset()
{
	// nothing to reset here
}

//-----------------------------------------------------------------------------
//! Calculates the stiffness matrix for 3D rigid elements.
//! Rigid elements don't generate stress, so there is nothing to do here
void FERigidSolidDomain::StiffnessMatrix(FELinearSystem& LS)
{
	// Caught you looking!
}

//-----------------------------------------------------------------------------
// Rigid bodies do not generate stress so there is nothing to do here
void FERigidSolidDomain::InternalForces(FEGlobalVector& R)
{
	// what you looking at ?!
}

//-----------------------------------------------------------------------------
//! We don't need to update the stresses for rigid elements
//!
void FERigidSolidDomain::Update(const FETimeInfo& tp)
{
	// Nothing to see here. Please move on.
}

//-----------------------------------------------------------------------------
void FERigidSolidDomain::MassMatrix(FELinearSystem& LS, double scale)
{
	// Only crickets here ... 
}

//-----------------------------------------------------------------------------
void FERigidSolidDomain::InertialForces(FEGlobalVector& R, std::vector<double>& F)
{
	// chirp, chirp ...
}

//-----------------------------------------------------------------------------
mat3d FERigidSolidDomain::CalculateMOI()
{
	mat3dd I(1); // identity tensor

	mat3d moi; moi.zero();

#pragma omp parallel 
	{
		vec3d r0[FEElement::MAX_NODES];
		mat3d moi_n; moi_n.zero();

		// loop over all elements
		int NE = Elements();
#pragma omp for nowait
		for (int iel = 0; iel < NE; ++iel)
		{
			FESolidElement& el = Element(iel);

			// initial coordinates
			int neln = el.Nodes();
			for (int i = 0; i < neln; ++i) r0[i] = GetMesh()->Node(el.m_node[i]).m_r0;

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
						moi_n += Iij * (H[i] * H[j] * Jw * dens);
					}
				}
			}
		}

#pragma omp critical
		moi += moi_n;
	}

	return moi;
}

double FERigidSolidDomain::CalculateMass()
{
	double mass = 0.0;
	int NE = Elements();
#pragma omp parallel for reduction(+:mass)
	for (int iel = 0; iel < NE; ++iel)
	{
		FESolidElement& el = Element(iel);

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

vec3d FERigidSolidDomain::CalculateCOM()
{
	vec3d rc(0, 0, 0);

	// loop over all elements
	int NE = Elements();
	#pragma omp parallel
	{
		vector<vec3d> r0(FEElement::MAX_NODES);
		vec3d rc_n(0, 0, 0);
#pragma omp for nowait
		for (int iel = 0; iel < NE; ++iel)
		{
			FESolidElement& el = Element(iel);

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

				// add to com and moi
				for (int i = 0; i < el.Nodes(); ++i)
				{
					rc_n += r0[i] * H[i] * detJ * gw[n] * dens;
				}
			}
		}

		#pragma omp critical
		rc += rc_n;
	}

	return rc;
}

//-----------------------------------------------------------------------------
void FERigidSolidDomain::BodyForce(FEGlobalVector& R, FEBodyForce& BF)
{
	// define some parameters that will be passed to lambda
	FEBodyForce* bodyForce = &BF;

	// evaluate the residual contribution
	LoadVector(R, m_dofU, [=](FEMaterialPoint& mp, int node_a, std::vector<double>& fa) {

		// evaluate density
		double density = m_pMat->Density(mp);

		// get the force
		vec3d f = bodyForce->force(mp);

		// get element shape functions
		double* H = mp.m_shape;

		// get the initial Jacobian
		double J0 = mp.m_J0;

		// set integrand
		fa[0] = -H[node_a] * density * f.x * J0;
		fa[1] = -H[node_a] * density * f.y * J0;
		fa[2] = -H[node_a] * density * f.z * J0;
	});
}
