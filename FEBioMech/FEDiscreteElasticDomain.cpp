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
#include "FEDiscreteElasticDomain.h"
#include "FEBioMech.h"
#include <FECore/FEMesh.h>
#include <FECore/FELinearSystem.h>

FEDiscreteElasticDomain::FEDiscreteElasticDomain(FEModel* fem) : FEDiscreteDomain(fem), FEElasticDomain(fem), m_dofU(fem), m_dofR(fem), m_dof(fem)
{
	m_pMat = nullptr;

	// TODO: Can this be done in Init, since there is no error checking
	if (fem)
	{
		m_dofU.AddVariable(FEBioMech::GetVariableName(FEBioMech::DISPLACEMENT));
		m_dofR.AddVariable(FEBioMech::GetVariableName(FEBioMech::RIGID_ROTATION));
	}
}

//-----------------------------------------------------------------------------
//! get the material (overridden from FEDomain)
FEMaterial* FEDiscreteElasticDomain::GetMaterial()
{
	return m_pMat;
}

//-----------------------------------------------------------------------------
//! set the material
void FEDiscreteElasticDomain::SetMaterial(FEMaterial* pmat)
{
	m_pMat = dynamic_cast<FEDiscreteElasticMaterial*>(pmat);
	assert(m_pMat);
}

//-----------------------------------------------------------------------------
// get the total dofs
const FEDofList& FEDiscreteElasticDomain::GetDOFList() const
{
	return m_dof;
}

//-----------------------------------------------------------------------------
void FEDiscreteElasticDomain::UnpackLM(FEElement &el, vector<int>& lm)
{
	int N = el.Nodes();
	lm.resize(N * 6);
	for (int i = 0; i<N; ++i)
	{
		FENode& node = m_pMesh->Node(el.m_node[i]);
		vector<int>& id = node.m_ID;

		// first the displacement dofs
		lm[3 * i] = id[m_dofU[0]];
		lm[3 * i + 1] = id[m_dofU[1]];
		lm[3 * i + 2] = id[m_dofU[2]];

		// rigid rotational dofs
		lm[3 * N + 3 * i] = id[m_dofR[0]];
		lm[3 * N + 3 * i + 1] = id[m_dofR[1]];
		lm[3 * N + 3 * i + 2] = id[m_dofR[2]];
	}
}

//-----------------------------------------------------------------------------
void FEDiscreteElasticDomain::Activate()
{
	for (int i = 0; i<Nodes(); ++i)
	{
		FENode& node = Node(i);
		if (node.HasFlags(FENode::EXCLUDE) == false)
		{
			if (node.m_rid < 0)
			{
				node.set_active(m_dofU[0]);
				node.set_active(m_dofU[1]);
				node.set_active(m_dofU[2]);
			}
		}
	}
}

//-----------------------------------------------------------------------------
void FEDiscreteElasticDomain::PreSolveUpdate(const FETimeInfo& timeInfo)
{
	FEDiscreteDomain::PreSolveUpdate(timeInfo);
	Update(timeInfo);
}

//-----------------------------------------------------------------------------
//! Calculates the forces due to discrete elements (i.e. springs)

void FEDiscreteElasticDomain::InternalForces(FEGlobalVector& R)
{
	FEMesh& mesh = *m_pMesh;

	vector<double> fe(6);
	vec3d u1, u2;

	vector<int> en(2), lm(6);

	for (size_t i = 0; i<m_Elem.size(); ++i)
	{
		// get the discrete element
		FEDiscreteElement& el = m_Elem[i];

		if (el.isActive())
		{
			// get the material point data
			FEMaterialPoint& mp = *el.GetMaterialPoint(0);
			FEDiscreteElasticMaterialPoint& ep = *mp.ExtractData<FEDiscreteElasticMaterialPoint>();

			// evaluate the force
			vec3d F = ep.m_Ft;

			// set up the force vector
			fe[0] = F.x;
			fe[1] = F.y;
			fe[2] = F.z;
			fe[3] = -F.x;
			fe[4] = -F.y;
			fe[5] = -F.z;

			// setup the node vector
			en[0] = el.m_node[0];
			en[1] = el.m_node[1];

			// set up the LM vector
			UnpackLM(el, lm);

			// assemble element
			R.Assemble(en, lm, fe);
		}
	}
}

//-----------------------------------------------------------------------------
//! Calculates the discrete element stiffness

void FEDiscreteElasticDomain::StiffnessMatrix(FELinearSystem& LS)
{
	FEMesh& mesh = *m_pMesh;

	FEElementMatrix ke;
	ke.resize(6, 6);
	ke.zero();

	vector<int> en(2), lm(6);

	// loop over all discrete elements
	for (size_t i = 0; i<m_Elem.size(); ++i)
	{
		// get the discrete element
		FEDiscreteElement& el = m_Elem[i];

		if (el.isActive())
		{
			// get the material point data
			FEDiscreteMaterialPoint& mp = *el.GetMaterialPoint(0)->ExtractData<FEDiscreteMaterialPoint>();

			// evaluate the stiffness
			mat3d A = m_pMat->Stiffness(mp);

			ke.zero();
			ke.add(0, 0, A); ke.add(0, 3, -A);
			ke.add(3, 0, -A); ke.add(3, 3, A);

			// setup the node vector
			en[0] = el.m_node[0];
			en[1] = el.m_node[1];

			// set up the LM vector
			UnpackLM(el, lm);

			// assemble the element into the global system
			ke.SetNodes(en);
			ke.SetIndices(lm);
			LS.Assemble(ke);
		}
	}
}

//-----------------------------------------------------------------------------
//! update domain data
void FEDiscreteElasticDomain::Update(const FETimeInfo& tp)
{
	FEMesh& mesh = *m_pMesh;

	double dt = tp.timeIncrement;

	// loop over all discrete elements
	for (size_t i = 0; i < m_Elem.size(); ++i)
	{
		// get the discrete element
		FEDiscreteElement& el = m_Elem[i];

		// only process active springs
		if (el.isActive())
		{
			// get the material point data
			FEDiscreteElasticMaterialPoint& mp = *el.GetMaterialPoint(0)->ExtractData<FEDiscreteElasticMaterialPoint>();

			// get the nodes of the element
			FENode& n1 = mesh.Node(el.m_node[0]);
			FENode& n2 = mesh.Node(el.m_node[1]);

			// get the nodal positions
			vec3d& rt1 = n1.m_rt;
			vec3d& rt2 = n2.m_rt;

			// get the previous nodal positions
			vec3d rp1 = n1.m_rp;
			vec3d rp2 = n2.m_rp;

			// get the initial nodal positions
			vec3d ri1 = n1.m_r0;
			vec3d ri2 = n2.m_r0;

			mp.m_drt = rt2 - rt1;
			mp.m_drp = rp2 - rp1;
			mp.m_dr0 = ri2 - ri1;

			mp.m_dvt = (mp.m_drt - mp.m_drp) / dt;

			// evaluate the force
			mp.m_Ft = m_pMat->Force(mp);
		}
	}
}
