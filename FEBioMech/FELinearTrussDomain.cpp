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
#include "FELinearTrussDomain.h"
#include <FECore/FEModel.h>
#include <FECore/FELinearSystem.h>
#include "FEBioMech.h"
#include <FECore/log.h>

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FELinearTrussDomain, FETrussDomain)
	ADD_PARAMETER(m_a0, "cross_sectional_area");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! Constructor
FELinearTrussDomain::FELinearTrussDomain(FEModel* pfem) : FETrussDomain(pfem), FEElasticDomain(pfem), m_dofU(pfem)
{
	m_pMat = 0;
	m_a0 = 0.0;
	m_dofU.AddVariable(FEBioMech::GetVariableName(FEBioMech::DISPLACEMENT));
}

//-----------------------------------------------------------------------------
//! copy operator
FELinearTrussDomain& FELinearTrussDomain::operator = (FELinearTrussDomain& d)
{ 
	m_Elem = d.m_Elem; 
	m_pMesh = d.m_pMesh; 
	return (*this); 
}

//-----------------------------------------------------------------------------
//! get the dof list
const FEDofList& FELinearTrussDomain::GetDOFList() const
{
	return m_dofU;
}

//-----------------------------------------------------------------------------
void FELinearTrussDomain::SetMaterial(FEMaterial* pmat)
{
	FETrussDomain::SetMaterial(pmat);
	m_pMat = dynamic_cast<FETrussMaterial*>(pmat);
	assert(m_pMat);
}

//-----------------------------------------------------------------------------
//! initialize the domain
bool FELinearTrussDomain::Init()
{
	if (m_a0 <= 0.0)
	{
		feLogError("Cross sectional area of \"%s\" must be positive.", GetName().c_str());
		return false;
	}

	for (int i = 0; i < Elements(); ++i)
	{
		FETrussElement& el = Element(i);
		el.m_a0 = m_a0;
	}
	
	return FETrussDomain::Init();
}

//-----------------------------------------------------------------------------
void FELinearTrussDomain::Reset()
{
	ForEachMaterialPoint([](FEMaterialPoint& mp) {
		mp.Init();
	});
}

//-----------------------------------------------------------------------------
void FELinearTrussDomain::UnpackLM(FEElement &el, vector<int>& lm)
{
	lm.resize(6);
	FENode& n1 = m_pMesh->Node(el.m_node[0]);
	FENode& n2 = m_pMesh->Node(el.m_node[1]);
	lm[0] = n1.m_ID[m_dofU[0]];
	lm[1] = n1.m_ID[m_dofU[1]];
	lm[2] = n1.m_ID[m_dofU[2]];
	lm[3] = n2.m_ID[m_dofU[0]];
	lm[4] = n2.m_ID[m_dofU[1]];
	lm[5] = n2.m_ID[m_dofU[2]];
}

//-----------------------------------------------------------------------------
void FELinearTrussDomain::Activate()
{
	for (int i=0; i<Nodes(); ++i)
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
void FELinearTrussDomain::PreSolveUpdate(const FETimeInfo& timeInfo)
{
	ForEachMaterialPoint([&](FEMaterialPoint& mp) {
		mp.Update(timeInfo);
	});
}

//-----------------------------------------------------------------------------

void FELinearTrussDomain::StiffnessMatrix(FELinearSystem& LS)
{
	int NT = (int)m_Elem.size();
	vector<int> lm;
	for (int iel =0; iel<NT; ++iel)
	{
		FETrussElement& el = m_Elem[iel];
		FEElementMatrix ke(el);
		ElementStiffness(iel, ke);
		UnpackLM(el, lm);
		ke.SetIndices(lm);
		LS.Assemble(ke);
	}
}

//-----------------------------------------------------------------------------
//! intertial stiffness matrix \todo implement this
void FELinearTrussDomain::MassMatrix(FELinearSystem& LS, double scale)
{
	for (int n = 0; n < Elements(); ++n)
	{
		FETrussElement& el = Element(n);

		matrix me(2, 2); me.zero();
		ElementMassMatrix(el, me);

		FEElementMatrix Me(6, 6); Me.zero();
		for (int i=0; i<2; ++i)
		{ 
			for (int j = 0; j < 2; ++j)
			{
				Me[3 * i    ][3 * j    ] = me[i][j];
				Me[3 * i + 1][3 * j + 1] = me[i][j];
				Me[3 * i + 2][3 * j + 2] = me[i][j];
			}
		}

		vector<int> lm;
		UnpackLM(el, lm);

		Me.SetIndices(lm);

		LS.Assemble(Me);
	}
}

//-----------------------------------------------------------------------------
void FELinearTrussDomain::ElementStiffness(int iel, matrix& ke)
{
	FETrussElement& el = Element(iel);

	// nodal coordinates
	vec3d r0[2], rt[2];
	for (int i=0; i<2; ++i)
	{
		r0[i] = m_pMesh->Node(el.m_node[i]).m_r0;
		rt[i] = m_pMesh->Node(el.m_node[i]).m_rt;
	}

	// intial length
	double L = (r0[1] - r0[0]).norm();

	// current length
	double l = (rt[1] - rt[0]).norm();

	// get the elastic tangent
	FEMaterialPoint& mp = *el.GetMaterialPoint(0);
	FETrussMaterialPoint& pt = *mp.ExtractData<FETrussMaterialPoint>();
	double E = m_pMat->Tangent(mp);

	// element initial volume
	double V = L*el.m_a0;

	// Kirchhoff Stress
	double tau = pt.m_tau;

	// scalar stiffness
	double k = V / (l*l)*( E - 2*tau);

	// axial force T = s*a = t*V/l
	double T = tau*V/l;

	// element normal
	vec3d n = TrussNormal(el);

	// calculate the tangent matrix
	ke.resize(6, 6);

	ke[0][0] = ke[3][3] = k*n.x*n.x + T/l;
	ke[1][1] = ke[4][4] = k*n.y*n.y + T/l;
	ke[2][2] = ke[5][5] = k*n.z*n.z + T/l;

	ke[0][1] = ke[1][0] = ke[3][4] = ke[4][3] = k*n.x*n.y;
	ke[1][2] = ke[2][1] = ke[4][5] = ke[5][4] = k*n.y*n.z;
	ke[0][2] = ke[2][0] = ke[3][5] = ke[5][3] = k*n.x*n.z;

	ke[0][3] = ke[3][0] = -ke[0][0]; ke[0][4] = ke[4][0] = -ke[0][1]; ke[0][5] = ke[5][0] = -ke[0][2];
	ke[1][3] = ke[3][1] = -ke[1][0]; ke[1][4] = ke[4][1] = -ke[1][1]; ke[1][5] = ke[5][1] = -ke[1][2];
	ke[2][3] = ke[3][2] = -ke[2][0]; ke[2][4] = ke[4][2] = -ke[2][1]; ke[2][5] = ke[5][2] = -ke[2][2];
}

//----------------------------------------------------------------------------
//! elemental mass matrix
void FELinearTrussDomain::ElementMassMatrix(FETrussElement& el, matrix& me)
{
	// nodal coordinates
	vec3d r0[2];
	for (int i = 0; i < 2; ++i)
	{
		r0[i] = m_pMesh->Node(el.m_node[i]).m_r0;
	}

	double L = (r0[1] - r0[0]).norm();
	double A = el.m_a0;
	double rho = m_pMat->Density();

	me[0][0] = rho * A * L / 3.0; me[0][1] = rho * A * L / 6.0;
	me[1][0] = rho * A * L / 6.0; me[1][1] = rho * A * L / 3.0;
}

//----------------------------------------------------------------------------
void FELinearTrussDomain::InternalForces(FEGlobalVector& R)
{
	// element force vector
	vector<double> fe;
	vector<int> lm;
	int NT = (int)m_Elem.size();
	for (int i=0; i<NT; ++i)
	{
		FETrussElement& el = m_Elem[i];
		ElementInternalForces(el, fe);
		UnpackLM(el, lm);
		R.Assemble(el.m_node, lm, fe);
	}
}

//-----------------------------------------------------------------------------
void FELinearTrussDomain::ElementInternalForces(FETrussElement& el, vector<double>& fe)
{
	FEMaterialPoint& mp = *el.GetMaterialPoint(0);
	FETrussMaterialPoint& pt = *(mp.ExtractData<FETrussMaterialPoint>());

	// get the element's normal
	vec3d n = TrussNormal(el);

	// get the element's Kirchhoff stress
	double tau = pt.m_tau;

	// nodal coordinates
	vec3d r0[2], rt[2];
	for (int i=0; i<2; ++i)
	{
		r0[i] = m_pMesh->Node(el.m_node[i]).m_r0;
		rt[i] = m_pMesh->Node(el.m_node[i]).m_rt;
	}

	// initial length
	double L = (r0[1] - r0[0]).norm();

	// current length
	double l = (rt[1] - rt[0]).norm();

	// elements initial volume
	double V = L*el.m_a0;

	// calculate nodal forces
	fe.resize(6);
	fe[0] = tau*V/l*n.x;
	fe[1] = tau*V/l*n.y;
	fe[2] = tau*V/l*n.z;
	fe[3] = -fe[0];
	fe[4] = -fe[1];
	fe[5] = -fe[2];
}

//-----------------------------------------------------------------------------
//! Update the truss' stresses
void FELinearTrussDomain::Update(const FETimeInfo& tp)
{
	// loop over all elements
	vec3d r0[2], rt[2];
	for (int i=0; i<(int) m_Elem.size(); ++i)
	{
		// unpack the element
		FETrussElement& el = m_Elem[i];

		// setup the material point
		FEMaterialPoint& mp = *(el.GetMaterialPoint(0));
		FETrussMaterialPoint& pt = *(mp.ExtractData<FETrussMaterialPoint>());

		// nodal coordinates
		for (int j=0; j<2; ++j)
		{
			r0[j] = m_pMesh->Node(el.m_node[j]).m_r0;
			rt[j] = m_pMesh->Node(el.m_node[j]).m_rt;
		}

		double l = (rt[1] - rt[0]).norm();
		double L = (r0[1] - r0[0]).norm();

		// calculate strain
		pt.m_l = l / L;

		// calculate stress
		pt.m_tau = m_pMat->Stress(mp);
	}
}
