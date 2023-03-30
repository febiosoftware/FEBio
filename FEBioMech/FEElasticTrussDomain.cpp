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
#include "FEElasticTrussDomain.h"
#include <FECore/FEModel.h>
#include <FECore/FELinearSystem.h>
#include "FEUncoupledMaterial.h"
#include "FEElasticMaterialPoint.h"
#include "FEBioMech.h"

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEElasticTrussDomain, FETrussDomain)
	ADD_PARAMETER(m_a0, "cross_sectional_area");
	ADD_PARAMETER(m_v, "v")->setLongName("Poisson's ratio");
END_FECORE_CLASS();
//-----------------------------------------------------------------------------
//! Constructor
FEElasticTrussDomain::FEElasticTrussDomain(FEModel* pfem) : FETrussDomain(pfem), FEElasticDomain(pfem), m_dofU(pfem)
{
	m_a0 = 0.0;
	m_v = 0.5;

	m_pMat = 0;
	// TODO: Can this be done in Init, since there is no error checking
	if (pfem)
	{
		m_dofU.AddVariable(FEBioMech::GetVariableName(FEBioMech::DISPLACEMENT));
	}
}

//-----------------------------------------------------------------------------
//! copy operator
FEElasticTrussDomain& FEElasticTrussDomain::operator = (FEElasticTrussDomain& d)
{ 
	m_Elem = d.m_Elem; 
	m_pMesh = d.m_pMesh; 
	return (*this); 
}

//-----------------------------------------------------------------------------
//! get the dof list
const FEDofList& FEElasticTrussDomain::GetDOFList() const
{
	return m_dofU;
}

//-----------------------------------------------------------------------------
void FEElasticTrussDomain::SetMaterial(FEMaterial* pmat)
{
	FETrussDomain::SetMaterial(pmat);
	m_pMat = dynamic_cast<FESolidMaterial*>(pmat);
	assert(m_pMat);
}

//-----------------------------------------------------------------------------
//! initialize the domain
bool FEElasticTrussDomain::Init()
{
	if (m_a0 != 0.0)
	{
		for (int i = 0; i < Elements(); ++i)
		{
			FETrussElement& el = Element(i);
			el.m_a0 = m_a0;
		}
	}


	for (int i = 0; i < (int)m_Elem.size(); ++i)
	{
		// unpack the element
		FETrussElement& el = m_Elem[i];

		// nodal coordinates
		vec3d r0[2];
		for (int j = 0; j < 2; ++j) r0[j] = m_pMesh->Node(el.m_node[j]).m_r0;

		// initial length
		el.m_L0 = (r0[1] - r0[0]).norm();
	}

	return FETrussDomain::Init();
}

//-----------------------------------------------------------------------------
void FEElasticTrussDomain::Reset()
{
	ForEachMaterialPoint([](FEMaterialPoint& mp) {
		mp.Init();
	});
}

//-----------------------------------------------------------------------------
void FEElasticTrussDomain::UnpackLM(FEElement &el, vector<int>& lm)
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
void FEElasticTrussDomain::Activate()
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
void FEElasticTrussDomain::PreSolveUpdate(const FETimeInfo& timeInfo)
{
	ForEachMaterialPoint([&](FEMaterialPoint& mp) {
		mp.Update(timeInfo);
	});
}

//-----------------------------------------------------------------------------

void FEElasticTrussDomain::StiffnessMatrix(FELinearSystem& LS)
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
//! intertial stiffness matrix
void FEElasticTrussDomain::MassMatrix(FELinearSystem& LS, double scale)
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
void FEElasticTrussDomain::ElementStiffness(int iel, matrix& ke)
{
	FETrussElement& el = Element(iel);

	// get the elastic tangent
	FEMaterialPoint& mp = *el.GetMaterialPoint(0);
	tens4ds c = m_pMat->Tangent(mp);
	double E = c(0, 0);

	// element initial volume
	double V = el.m_L0*el.m_a0;
	double l = el.m_L0 / el.m_lam;

	// Kirchhoff Stress
	double tau = el.m_tau;

	// scalar stiffness
	double k = E;// V / (l * l) * (E - 2 * tau);

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
void FEElasticTrussDomain::ElementMassMatrix(FETrussElement& el, matrix& me)
{
	// nodal coordinates
	vec3d r0[2];
	for (int i = 0; i < 2; ++i)
	{
		r0[i] = m_pMesh->Node(el.m_node[i]).m_r0;
	}

	double L = (r0[1] - r0[0]).norm();
	double A = el.m_a0;
	double rho = m_pMat->Density(*el.GetMaterialPoint(0));

	me[0][0] = rho * A * L / 3.0; me[0][1] = rho * A * L / 6.0;
	me[1][0] = rho * A * L / 6.0; me[1][1] = rho * A * L / 3.0;
}

//----------------------------------------------------------------------------
void FEElasticTrussDomain::InternalForces(FEGlobalVector& R)
{
	// element force vector
	int NT = (int)m_Elem.size();
#pragma omp parallel for
	for (int i=0; i<NT; ++i)
	{
		FETrussElement& el = m_Elem[i];

		vector<double> fe;
		ElementInternalForces(el, fe);

		vector<int> lm;
		UnpackLM(el, lm);
		R.Assemble(el.m_node, lm, fe);
	}
}

//-----------------------------------------------------------------------------
void FEElasticTrussDomain::ElementInternalForces(FETrussElement& el, vector<double>& fe)
{
	FEMaterialPoint& mp = *el.GetMaterialPoint(0);
	FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());

	// get the element's normal
	vec3d n = TrussNormal(el);

	// nodal coordinates
	vec3d r0[2] = {
		m_pMesh->Node(el.m_node[0]).m_r0,
		m_pMesh->Node(el.m_node[1]).m_r0
	};

	vec3d rt[2] = {
		m_pMesh->Node(el.m_node[0]).m_rt,
		m_pMesh->Node(el.m_node[1]).m_rt
	};

	// initial length
	double L = (r0[1] - r0[0]).norm();

	// current length
	double l = (rt[1] - rt[0]).norm();

	// stress in element
	double tau = pt.m_s.xx();

	// elements initial volume
	double V = L*el.m_a0;

	// current volume
	double v = V * pt.m_J;

	// current area
	double a = v / l;

	// calculate nodal forces
	fe.resize(6);
	fe[0] = tau*a*n.x;
	fe[1] = tau*a*n.y;
	fe[2] = tau*a*n.z;
	fe[3] = -fe[0];
	fe[4] = -fe[1];
	fe[5] = -fe[2];
}

//-----------------------------------------------------------------------------
//! Update the truss' stresses
void FEElasticTrussDomain::Update(const FETimeInfo& tp)
{
	FEUncoupledMaterial* um = dynamic_cast<FEUncoupledMaterial*>(GetMaterial());

	// loop over all elements
#pragma omp parallel for
	for (int i=0; i<(int) m_Elem.size(); ++i)
	{
		// unpack the element
		FETrussElement& el = m_Elem[i];

		// nodal coordinates
		vec3d r0[2], rt[2];
		for (int j=0; j<2; ++j)
		{
			r0[j] = m_pMesh->Node(el.m_node[j]).m_r0;
			rt[j] = m_pMesh->Node(el.m_node[j]).m_rt;
		}

		double l = (rt[1] - rt[0]).norm();
		double L = (r0[1] - r0[0]).norm();

		double lam = l / L;
		double linv2 = 1.0 / sqrt(pow(lam, 2.0*m_v));

		// calculate strain
		el.m_lam = l / L;

		// setup deformation gradient (assuming incompressibility)
		mat3d F(lam, 0.0, 0.0, 0.0, linv2, 0.0, 0.0, 0.0, linv2);

		// calculate stress
		FEMaterialPoint& mp = *el.GetMaterialPoint(0);
		FEElasticMaterialPoint& ep = *mp.ExtractData<FEElasticMaterialPoint>();
		ep.m_F = F;
		ep.m_J = F.det();

		if (um)
		{
			mat3ds devs = um->DevStress(mp);
			double p = -devs.yy();
			ep.m_s = devs + mat3dd(p);
		}
		else
		{
			ep.m_s = m_pMat->Stress(mp);
		}
	}
}
