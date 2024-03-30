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
#include "FEBodyForce.h"
#include "FEBioMech.h"
#include <FECore/log.h>

BEGIN_FECORE_CLASS(FELinearTrussDomain, FETrussDomain)
	ADD_PARAMETER(m_a0, "cross_sectional_area");
END_FECORE_CLASS();

FELinearTrussDomain::FELinearTrussDomain(FEModel* pfem) : FETrussDomain(pfem), FEElasticDomain(pfem), m_dofU(pfem), m_dofV(pfem), m_dofR(pfem)
{
	m_pMat = nullptr;
	m_a0 = 0.0;
	m_dofU.AddVariable(FEBioMech::GetVariableName(FEBioMech::DISPLACEMENT));
	m_dofR.AddVariable(FEBioMech::GetVariableName(FEBioMech::RIGID_ROTATION));
	m_dofV.AddVariable(FEBioMech::GetVariableName(FEBioMech::VELOCITY));
}

const FEDofList& FELinearTrussDomain::GetDOFList() const
{
	return m_dofU;
}

void FELinearTrussDomain::SetMaterial(FEMaterial* pmat)
{
	FETrussDomain::SetMaterial(pmat);
	m_pMat = dynamic_cast<FETrussMaterial*>(pmat);
	assert(m_pMat);
}

bool FELinearTrussDomain::Init()
{
	if (m_a0 <= 0.0)
	{
		feLogError("Cross sectional area of \"%s\" must be positive.", GetName().c_str());
		return false;
	}

	for (FETrussElement& el : m_Elem)
	{
		el.m_a0 = m_a0;
	}
	
	return FETrussDomain::Init();
}

void FELinearTrussDomain::Reset()
{
	ForEachMaterialPoint([](FEMaterialPoint& mp) {
		mp.Init();
	});
}

void FELinearTrussDomain::UnpackLM(FEElement &el, vector<int>& lm)
{
	lm.resize(12);
	FENode& n1 = m_pMesh->Node(el.m_node[0]);
	FENode& n2 = m_pMesh->Node(el.m_node[1]);
	lm[0] = n1.m_ID[m_dofU[0]];
	lm[1] = n1.m_ID[m_dofU[1]];
	lm[2] = n1.m_ID[m_dofU[2]];
	lm[3] = n2.m_ID[m_dofU[0]];
	lm[4] = n2.m_ID[m_dofU[1]];
	lm[5] = n2.m_ID[m_dofU[2]];

	// rigid rotational dofs
	lm[6 ] = n1.m_ID[m_dofR[0]];
	lm[7 ] = n1.m_ID[m_dofR[1]];
	lm[8 ] = n1.m_ID[m_dofR[2]];
	lm[9 ] = n2.m_ID[m_dofR[0]];
	lm[10] = n2.m_ID[m_dofR[1]];
	lm[11] = n2.m_ID[m_dofR[2]];
}

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

void FELinearTrussDomain::PreSolveUpdate(const FETimeInfo& timeInfo)
{
	ForEachMaterialPoint([&](FEMaterialPoint& mp) {
		mp.Update(timeInfo);
	});
}

void FELinearTrussDomain::MassMatrix(FELinearSystem& LS, double scale)
{
	for (FETrussElement& el : m_Elem)
	{
		matrix me(2, 2); me.zero();
		ElementMassMatrix(el, me);

		// The element mass matrix calculated above only evaluates the 2x2 mass matrix
		// Thus, we need to inlate it to a 6x6 mass matrix.
		FEElementMatrix Me(6, 6); Me.zero();
		for (int i=0; i<2; ++i)
		{ 
			for (int j = 0; j < 2; ++j)
			{
				Me[3 * i    ][3 * j    ] = scale*me[i][j];
				Me[3 * i + 1][3 * j + 1] = scale*me[i][j];
				Me[3 * i + 2][3 * j + 2] = scale*me[i][j];
			}
		}

		vector<int> lm;
		UnpackLM(el, lm);
		Me.SetNodes(el.m_node);
		Me.SetIndices(lm);
		LS.Assemble(Me);
	}
}

void FELinearTrussDomain::ElementMassMatrix(FETrussElement& el, matrix& me)
{
	double L = el.m_L0;
	double A = el.m_a0;
	double rho = m_pMat->Density();

	me[0][0] = rho * A * L / 3.0; me[0][1] = rho * A * L / 6.0;
	me[1][0] = rho * A * L / 6.0; me[1][1] = rho * A * L / 3.0;
}

void FELinearTrussDomain::StiffnessMatrix(FELinearSystem& LS)
{
	vector<int> lm;
	for (FETrussElement& el : m_Elem)
	{
		FEElementMatrix ke(el);
		ElementStiffness(el, ke);
		UnpackLM(el, lm);
		ke.SetNodes(el.m_node);
		ke.SetIndices(lm);
		LS.Assemble(ke);
	}
}

void FELinearTrussDomain::ElementStiffness(FETrussElement& el, matrix& ke)
{
	// intial length
	double L = el.m_L0;

	// current length
	double l = el.m_Lt;

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

	vec3d n = GetTrussAxisVector(el);

	// calculate the tangent matrix
	ke.resize(6, 6);

	FETimeInfo& ti = GetFEModel()->GetTime();
	double alphaf = ti.alphaf;

	ke[0][0] = ke[3][3] = (k*n.x*n.x + T/l)*alphaf;
	ke[1][1] = ke[4][4] = (k*n.y*n.y + T/l)*alphaf;
	ke[2][2] = ke[5][5] = (k*n.z*n.z + T/l)*alphaf;

	ke[0][1] = ke[1][0] = ke[3][4] = ke[4][3] = k*n.x*n.y*alphaf;
	ke[1][2] = ke[2][1] = ke[4][5] = ke[5][4] = k*n.y*n.z*alphaf;
	ke[0][2] = ke[2][0] = ke[3][5] = ke[5][3] = k*n.x*n.z*alphaf;

	ke[0][3] = ke[3][0] = -ke[0][0]; ke[0][4] = ke[4][0] = -ke[0][1]; ke[0][5] = ke[5][0] = -ke[0][2];
	ke[1][3] = ke[3][1] = -ke[1][0]; ke[1][4] = ke[4][1] = -ke[1][1]; ke[1][5] = ke[5][1] = -ke[1][2];
	ke[2][3] = ke[3][2] = -ke[2][0]; ke[2][4] = ke[4][2] = -ke[2][1]; ke[2][5] = ke[5][2] = -ke[2][2];
}

void FELinearTrussDomain::InternalForces(FEGlobalVector& R)
{
	// element force vector
	vector<double> fe;
	vector<int> lm;
	for (FETrussElement& el : m_Elem)
	{
		ElementInternalForces(el, fe);
		UnpackLM(el, lm);
		R.Assemble(el.m_node, lm, fe);
	}
}

void FELinearTrussDomain::ElementInternalForces(FETrussElement& el, vector<double>& fe)
{
	FEMaterialPoint& mp = *el.GetMaterialPoint(0);
	FETrussMaterialPoint& pt = *(mp.ExtractData<FETrussMaterialPoint>());

	FETimeInfo& ti = GetFEModel()->GetTime();
	double alphaf = ti.alphaf;

	vec3d rt[2];
	for (int i = 0; i < 2; ++i)
	{
		FENode& nd = m_pMesh->Node(el.m_node[i]);
		rt[i] = nd.m_rt * alphaf + nd.m_rp * (1 - alphaf);
	}
	vec3d n = rt[1] - rt[0];
	n.unit();

	// get the element's Kirchhoff stress
	double tau = pt.m_tau;

	// initial length
	double L = el.m_L0;

	// current length
	double l = el.m_Lt;

	// initial volume
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

//! Calculates inertial forces for dynamic analyses
void FELinearTrussDomain::InertialForces(FEGlobalVector& R, vector<double>& F)
{
	vector<double> fe;
	vector<int> lm;
	for (FETrussElement& el : m_Elem)
	{
		ElementInertialForces(el, fe);
		UnpackLM(el, lm);
		R.Assemble(el.m_node, lm, fe);
	}
}

void FELinearTrussDomain::ElementInertialForces(FETrussElement& el, vector<double>& fe)
{
	FEMaterialPoint& mp = *el.GetMaterialPoint(0);
	FETrussMaterialPoint& pt = *(mp.ExtractData<FETrussMaterialPoint>());

	FETimeInfo& ti = GetFEModel()->GetTime();
	double alpham = ti.alpham;

	// get the nodal accelerations
	vec3d a[2];
	int neln = el.Nodes();
	for (int i = 0; i < neln; ++i)
	{
		FENode& node = m_pMesh->Node(el.m_node[i]);
		a[i] = node.m_ap* (1 - alpham) + node.m_at* alpham;
	}

	// calculate element contribution to inertial force
	matrix me(2, 2);
	ElementMassMatrix(el, me);

	fe.resize(6, 0.0);
	fe[0] = -(me[0][0] * a[0].x + me[0][1] * a[1].x);
	fe[1] = -(me[0][0] * a[0].y + me[0][1] * a[1].y);
	fe[2] = -(me[0][0] * a[0].z + me[0][1] * a[1].z);
	fe[3] = -(me[1][0] * a[0].x + me[1][1] * a[1].x);
	fe[4] = -(me[1][0] * a[0].y + me[1][1] * a[1].y);
	fe[5] = -(me[1][0] * a[0].z + me[1][1] * a[1].z);
}

//! calculate body force on a truss domain
void FELinearTrussDomain::BodyForce(FEGlobalVector& R, FEBodyForce& bf)
{
	vector<double> fe(6);
	vector<int> lm;

	double density = m_pMat->Density();

	for (FETrussElement& el : m_Elem)
	{
		zero(fe);
		int nint = el.GaussPoints();
		double* w = el.GaussWeights();
		for (int n = 0; n < nint; ++n)
		{
			FEMaterialPoint& mp = *el.GetMaterialPoint(n);

			// get the force
			vec3d f = bf.force(mp);

			// get element shape functions
			double* H = el.H(n);

			// get the initial Jacobian
			double J0 = el.m_a0 * el.m_L0 / 2.0;

			// set integrand
			fe[0] -= H[0] * density * f.x * J0 * w[n];
			fe[1] -= H[0] * density * f.y * J0 * w[n];
			fe[2] -= H[0] * density * f.z * J0 * w[n];
			fe[3] -= H[1] * density * f.x * J0 * w[n];
			fe[4] -= H[1] * density * f.y * J0 * w[n];
			fe[5] -= H[1] * density * f.z * J0 * w[n];
		}

		UnpackLM(el, lm);
		R.Assemble(el.m_node, lm, fe);
	}
}

//! body force stiffness matrix
void FELinearTrussDomain::BodyForceStiffness(FELinearSystem& LS, FEBodyForce& bf)
{
	double density = m_pMat->Density();

	vector<int> lm;
	for (FETrussElement& el : m_Elem)
	{
		FEElementMatrix ke(6, 6); ke.zero();
		int nint = el.GaussPoints();
		int neln = el.Nodes();
		double* gw = el.GaussWeights();
		for (int n = 0; n < nint; ++n)
		{
			FEMaterialPoint& mp = *el.GetMaterialPoint(n);

			double* H = el.H(n);

			double J0 = el.m_a0 * el.m_L0 / 2.0;

			double Jdw = J0 * density * gw[n];

			mat3d K = bf.stiffness(mp);

			// put it together
			for (int a=0; a<neln; ++a)
				for (int b = 0; b < neln; ++b)
				{
					ke.sub(3 * a, 3 * b, K * (H[a] * H[b] * Jdw));
				}
		}

		UnpackLM(el, lm);
		ke.SetIndices(lm);
		ke.SetNodes(el.m_node);
		LS.Assemble(ke);
	}
}

void FELinearTrussDomain::Update(const FETimeInfo& tp)
{
	// get the poisson's ratio. 
	FELinearTrussMaterial* mat = dynamic_cast<FELinearTrussMaterial*>(m_pMat);
	double nu = (mat ? mat->m_v : 0.5);

	double alpham = tp.alpham;
	double alphaf = tp.alphaf;

	FEMesh& mesh = *m_pMesh;
	for (FETrussElement& el : m_Elem)
	{
		// calculate the current length
		vec3d rt[2], vt[2], at[2];
		for (int j = 0; j < el.Nodes(); ++j)
		{
			FENode& node = m_pMesh->Node(el.m_node[j]);
			rt[j] = node.m_rt * alphaf + node.m_rp*(1 - alphaf);
			vec3d vtj = node.get_vec3d(m_dofV[0], m_dofV[1], m_dofV[2]);
			vt[j] = vtj * alphaf + node.m_vp * (1 - alphaf);
			at[j] = node.m_at * alpham + node.m_ap * (1 - alpham);
		}

		el.m_Lt = (rt[1] - rt[0]).norm();
		vec3d e = (rt[1] - rt[0]).normalized();

		for (int n = 0; n < el.GaussPoints(); ++n)
		{
			FEMaterialPoint& mp = *(el.GetMaterialPoint(n));
			FETrussMaterialPoint& pt = *(mp.ExtractData<FETrussMaterialPoint>());

			mp.m_rt = el.Evaluate(rt, n);
			pt.m_v = el.Evaluate(vt, n);
			pt.m_a = el.Evaluate(at, n);

			// calculate stretch
			pt.m_lam = el.m_Lt / el.m_L0;

			// volume ratio (this assume linear truss material!)
			double J = pow(pt.m_lam, 1.0 - 2.0 * nu);

			// calculate stress
			pt.m_tau = m_pMat->Stress(mp);
			pt.m_s = dyad(e) * (pt.m_tau / J);
		}
	}
}
