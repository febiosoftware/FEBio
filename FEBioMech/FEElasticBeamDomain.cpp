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
#include "FEBodyForce.h"
#include <FECore/FEMesh.h>
#include <FECore/FEModel.h>
#include <FECore/FEAnalysis.h>
#include <FECore/FESolver.h>
#include <FECore/fecore_debug.h>


FEElasticBeamDomain::FEElasticBeamDomain(FEModel* fem) : FEBeamDomain(fem), FEElasticDomain(fem), 
	m_dofs(fem), m_dofQ(fem), m_dofV(fem), m_dofW(fem), m_dofA(fem)
{
	m_mat = nullptr;

	if (fem)
	{
		m_dofs.AddVariable(FEBioMech::GetVariableName(FEBioMech::DISPLACEMENT));
		m_dofs.AddVariable(FEBioMech::GetVariableName(FEBioMech::ROTATION));

		m_dofQ.AddVariable(FEBioMech::GetVariableName(FEBioMech::ROTATION));
		m_dofV.AddVariable(FEBioMech::GetVariableName(FEBioMech::VELOCITY));
		m_dofW.AddVariable(FEBioMech::GetVariableName(FEBioMech::BEAM_ANGULAR_VELOCITY));
		m_dofA.AddVariable(FEBioMech::GetVariableName(FEBioMech::BEAM_ANGULAR_ACCELERATION));
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
	int etype = (espec.eshape == ET_LINE3 ? FE_BEAM3G2 : FE_BEAM2G2);
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

		// construct beam coordinate system
		vec3d E3 = r0[1] - r0[0]; E3.Normalize();
		vec3d E1, E2;
		if (fabs(E3 * vec3d(1, 0, 0)) > 0.5) E2 = E3 ^ vec3d(0, 1, 0);
		else E2 = E3 ^ vec3d(1,0,0);
		E2.Normalize();
		E1 = E2 ^ E3;
		el.m_E = mat3d(E1, E2, E3);
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

		vector<int> lm(6*ne);
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

		// copy local coordinate system (Probably don't need to do this every time)
		mp.m_Q = el.m_E;

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
	// We need this for the incremental update of prescribed rotational dofs
	int niter = GetFEModel()->GetCurrentStep()->GetFESolver()->m_niter;

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
		// NOTE: The ui vector does not contain the prescribed displacement updates!!
		vector<vec3d> dri(ne);
		for (int j = 0; j < ne; ++j)
		{
			FENode& node = Node(el.m_lnode[j]);
			std::vector<int>& id = node.m_ID;
			int eq[3] = { id[m_dofs[3]], id[m_dofs[4]], id[m_dofs[5]] };
			vec3d d(0, 0, 0);
			if (eq[0] >= 0) d.x = ui[eq[0]];
			if (eq[1] >= 0) d.y = ui[eq[1]];
			if (eq[2] >= 0) d.z = ui[eq[2]];

			// for prescribed nodes, determine the rotation increment
			if ((eq[0] < 0) && (eq[1] < 0) && (eq[2] < 0))
			{
				if (niter == 0)
				{
					vec3d Rp, Rt;
					int n;
					n = -eq[0] - 1; if (n >= 0) { Rt.x = node.get(m_dofs[3]); Rp.x = node.get_prev(m_dofs[3]); }
					n = -eq[1] - 1; if (n >= 0) { Rt.y = node.get(m_dofs[4]); Rp.y = node.get_prev(m_dofs[4]); }
					n = -eq[2] - 1; if (n >= 0) { Rt.z = node.get(m_dofs[5]); Rp.z = node.get_prev(m_dofs[5]); }
					quatd Qp(Rp), Qt(Rt);
					quatd dQ = Qp.Conjugate() * Qt;
					d = dQ.GetRotationVector();
				}
			}

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
			mat3da Wn(mp.m_kn);
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
			mp.m_k = w1 + w2;

			if (finalFlag)
			{
				mp.m_Ri = dq * mp.m_Ri;
				mp.m_kn = mp.m_k;
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
	constexpr int NMAX = FEElement::MAX_NODES;
	vec3d rt[NMAX], vt[NMAX], at[NMAX], wt[NMAX], alt[NMAX];
	int ne = el.Nodes();
	for (int i = 0; i < ne; ++i)
	{
		FENode& node = Node(el.m_lnode[i]);
		rt[i] = node.m_rt;
		vt[i] = node.get_vec3d(m_dofV[0], m_dofV[1], m_dofV[2]);
		at[i] = node.m_at;

		wt[i]  = node.get_vec3d(m_dofW[0], m_dofW[1], m_dofW[2]);
		alt[i] = node.get_vec3d(m_dofA[0], m_dofA[1], m_dofA[2]);
	}

	// initial length
	double L0 = el.m_L0;
	double J = L0 / 2;

	FEElasticBeamMaterial& mat = *m_mat;
	double rho = mat.m_density;
	double A_rho = mat.m_A * rho;

	double I1 = rho * mat.m_I1;
	double I2 = rho * mat.m_I2;
	double I3 = I1 + I2;

	vec3d E1 = el.m_E.col(0);
	vec3d E2 = el.m_E.col(1);
	vec3d E3 = el.m_E.col(2);

	mat3d I0 = (E1 & E1)*I1 + (E2 & E2)*I2 + (E3 & E3) * I3; // material inertia tensor

	// loop over all integration points
	int nint = el.GaussPoints();
	for (int n = 0; n < nint; ++n)
	{
		// get the material point
		FEElasticBeamMaterialPoint& mp = *(el.GetMaterialPoint(n)->ExtractData<FEElasticBeamMaterialPoint>());

		// update quantities for dynamics
		mp.m_vt = el.Evaluate(vt, n);
		mp.m_at = el.Evaluate(at, n);
		mp.m_dpt = mp.m_at * A_rho;

		// (spatial) rotational quantities
		vec3d w = el.Evaluate(wt, n);
		vec3d al = el.Evaluate(alt, n);
		mp.m_wt = w;
		mp.m_alt = al;

		// calculate the spatial inertia tensor
		mat3d R = mp.m_Rt.RotationMatrix();
		mat3d I = R * I0 * R.transpose();

		// calculate rate of angular momentum
		mp.m_dht = I * al + (w ^ (I * w));

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
		mp.m_Gamma = qi * G0 - E3;
		mp.m_Kappa = qi * mp.m_k; // m_k is updated in IncrementalUpdate(std::vector<double>& ui)

		// evaluate the (spatial) stress
		mp.m_Q = el.m_E;
		m_mat->Stress(mp);
	}
}

//! Calculates the inertial forces (for dynamic problems)
void FEElasticBeamDomain::InertialForces(FEGlobalVector& R, std::vector<double>& F)
{
	for (auto& el : m_Elem)
	{
		int neln = el.Nodes();
		vector<double> fe(6 * neln, 0.0);
		ElementInertialForce(el, fe);
		vector<int> lm(6 * neln);
		UnpackLM(el, lm);
		R.Assemble(lm, fe);
	}
}

void FEElasticBeamDomain::ElementInertialForce(FEBeamElement& el, std::vector<double>& fe)
{
	int nint = el.GaussPoints();
	int neln = el.Nodes();
	double* gw = el.GaussWeights();
	for (int n = 0; n < nint; ++n)
	{
		FEElasticBeamMaterialPoint& mp = *(el.GetMaterialPoint(n)->ExtractData<FEElasticBeamMaterialPoint>());

		double* H = el.H(n);

		double J = el.m_L0 / 2.0;
		double Jw = J * gw[n];

		vec3d dp = mp.m_dpt;
		vec3d dh = mp.m_dht;

		for (int i = 0; i < neln; ++i)
		{
			fe[6*i    ] -= H[i]*dp.x*Jw;
			fe[6*i + 1] -= H[i]*dp.y*Jw;
			fe[6*i + 2] -= H[i]*dp.z*Jw;
			fe[6*i + 3] -= H[i]*dh.x*Jw;
			fe[6*i + 4] -= H[i]*dh.y*Jw;
			fe[6*i + 5] -= H[i]*dh.z*Jw;
		}
	}
}

//! Calculates the mass matrix (for dynamic problems)
void FEElasticBeamDomain::MassMatrix(FELinearSystem& LS, double scale) 
{
	for (auto& el : m_Elem)
	{
		int neln = el.Nodes();
		FEElementMatrix ke(6 * neln, 6 * neln); ke.zero();
		ElementMassMatrix(el, ke);
		vector<int> lm(6 * neln);
		UnpackLM(el, lm);
		ke.SetIndices(lm);
		ke.SetNodes(el.m_node);
		LS.Assemble(ke);
	}
}

// tanx = tan(x)/x
double tanx(double x)
{
	double r = (fabs(x) < 1e-9 ? 1.0 : tan(x) / x);
	return r;
}

void FEElasticBeamDomain::ElementMassMatrix(FEBeamElement& el, FEElementMatrix& ke)
{
	int nint = el.GaussPoints();
	int neln = el.Nodes();
	double* gw = el.GaussWeights();

	FEElasticBeamMaterial& mat = *m_mat;
	double rho = mat.m_density;
	double A_rho = mat.m_A * rho;

	double I1 = rho * mat.m_I1;
	double I2 = rho * mat.m_I2;
	double I3 = I1 + I2;

	vec3d E1 = el.m_E.col(0);
	vec3d E2 = el.m_E.col(1);
	vec3d E3 = el.m_E.col(2);

	mat3d I0 = (E1 & E1) * I1 + (E2 & E2) * I2 + (E3 & E3) * I3; // material inertia tensor

	FETimeInfo& ti = GetFEModel()->GetTime();
	double h = ti.timeIncrement;
	double b = ti.beta;
	double g = ti.gamma;
	double h2bi = 1.0 / (h * h * b);
	double hg = h * g;

	// When evaluating at time 0 (to determine initial accelerations), we need 
	// to make a minor change. 
	if (ti.currentTime == 0.0) { h2bi = hg = 1.0; }

	for (int n = 0; n < nint; ++n)
	{
		FEElasticBeamMaterialPoint& mp = *(el.GetMaterialPoint(n)->ExtractData<FEElasticBeamMaterialPoint>());

		vec3d ri = mp.m_Ri.GetRotationVector();
		vec3d e = ri.Normalized();
		double th = ri.norm();
		mat3da ts(ri);
		mat3dd I(1.0);
		mat3ds exe = dyad(e);
		mat3d T = exe + (I - exe)*(1.0/tanx(th*0.5)) - ts*0.5;

		double* H = el.H(n);

		double J = el.m_L0 / 2.0;
		double Jw = J * gw[n];

		mat3d Rt = mp.m_Rt.RotationMatrix();
		mat3d Rp = mp.m_Rp.RotationMatrix();

		vec3d Wt = Rt.transpose() * mp.m_wt;
		vec3d At = Rt.transpose() * mp.m_alt;
		mat3da W_hat(Wt);
		mat3da Hs(I0 * At + (Wt ^ (I0 * Wt)));
		mat3da IW(I0 * Wt);

		double M11 = A_rho * Jw * h2bi;
		mat3d M22 = (-Rt * Hs + Rt * (I0 - IW * hg + (W_hat * I0) * hg)*h2bi)* Rp.transpose()* T;

		for (int i=0; i<neln; ++i)
			for (int j = 0; j < neln; ++j)
			{
				mat3d m1; m1.zero();
				m1[0][0] = M11 * H[i] * H[j];
				m1[1][1] = M11 * H[i] * H[j];
				m1[2][2] = M11 * H[i] * H[j];

				mat3d m2;
				m2 = M22 * (H[i] * H[j] * Jw);

				int I = 6*i, J = 6*j;
				ke.add(I, J, m1);
				ke.add(I + 3, J + 3, m2);
			}
	}
}

void FEElasticBeamDomain::BodyForce(FEGlobalVector& R, FEBodyForce& bf) 
{
	for (int iel = 0; iel < Elements(); ++iel)
	{
		FEBeamElement& el = Element(iel);

		int ne = el.Nodes();
		int ndof = ne * 6;
		vector<double> fe(ndof, 0.0);

		ElementBodyForce(el, fe, bf);

		vector<int> lm(6 * ne);
		UnpackLM(el, lm);
		R.Assemble(lm, fe);
	}
}

void FEElasticBeamDomain::ElementBodyForce(FEBeamElement& el, std::vector<double>& fe, FEBodyForce& bf)
{
	// reference length of beam
	double L0 = el.m_L0;

	FEElasticBeamMaterial& mat = *m_mat;
	double rho = mat.m_density;
	double A_rho = mat.m_A * rho;

	// loop over integration points
	int nint = el.GaussPoints();
	int ne = el.Nodes();
	double* w = el.GaussWeights();
	for (int n = 0; n < nint; ++n)
	{
		FEMaterialPoint& mp = *el.GetMaterialPoint(n);
		FEElasticBeamMaterialPoint& ebm = *(mp.ExtractData<FEElasticBeamMaterialPoint>());

		// shape functions
		double* H = el.H(n);

		// J = dS / dr
		double J = L0 / 2.0;

		double wJA = w[n] * J * A_rho;

		vec3d fn = bf.force(mp);

		for (int i = 0; i < ne; ++i)
		{
			fe[i * 6    ] += H[i] * fn.x * wJA;
			fe[i * 6 + 1] += H[i] * fn.y * wJA;
			fe[i * 6 + 2] += H[i] * fn.z * wJA;
		}
	}
}

void FEElasticBeamDomain::BodyForceStiffness(FELinearSystem& LS, FEBodyForce& bf) 
{
	for (int iel = 0; iel < Elements(); ++iel)
	{
		FEBeamElement& el = Element(iel);

		int ne = el.Nodes();
		int ndof = ne * 6;

		FEElementMatrix ke(ndof, ndof);
		ke.zero();
		ElementBodyForceStiffness(el, ke, bf);

		vector<int> lm(ndof);
		UnpackLM(el, lm);
		ke.SetIndices(lm);
		LS.Assemble(ke);
	}
}

void FEElasticBeamDomain::ElementBodyForceStiffness(FEBeamElement& el, FEElementMatrix& ke, FEBodyForce& bf)
{
	// reference length of beam
	double L0 = el.m_L0;

	FEElasticBeamMaterial& mat = *m_mat;
	double rho = mat.m_density;
	double A_rho = mat.m_A * rho;

	// loop over integration points
	int nint = el.GaussPoints();
	int ne = el.Nodes();
	double* w = el.GaussWeights();
	for (int n = 0; n < nint; ++n)
	{
		FEMaterialPoint& mp = *el.GetMaterialPoint(n);
		FEElasticBeamMaterialPoint& ebm = *(mp.ExtractData<FEElasticBeamMaterialPoint>());

		// shape functions
		double* H = el.H(n);

		// J = dS / dr
		double J = L0 / 2.0;

		double wJA = w[n] * J * A_rho;

		mat3d k = bf.stiffness(mp);

		for (int a = 0; a < ne; ++a)
			for (int b = 0; b < ne; ++b)
			{
				for (int i=0; i<3; ++i)
					for (int j = 0; j < 3; ++j)
					{
						ke(a * 6 + i, b * 6 + j) -= H[a] * k[i][j] * H[b] * wJA;
					}
			}
	}
}
