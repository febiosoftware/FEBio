/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2020 University of Utah, The Trustees of Columbia University in
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
#include "FEContinuousElasticDamage.h"

//=========================================================================================

class FEContinuousElasticDamage::Data : public FEMaterialPoint
{
public:
	Data(FEMaterialPoint* pm) : FEMaterialPoint(pm) {}
};

//=========================================================================================

BEGIN_FECORE_CLASS(FEContinuousElasticDamage, FEElasticMaterial)
	ADD_PARAMETER(m_c1, FE_RANGE_GREATER(0.0), "c1");
	ADD_PARAMETER(m_k, FE_RANGE_GREATER(0.0), "k");
	ADD_PARAMETER(m_a1, FE_RANGE_GREATER_OR_EQUAL(0.0), "a1");
	ADD_PARAMETER(m_a2, FE_RANGE_GREATER(1.0), "a2");
	ADD_PARAMETER(m_kappa, FE_RANGE_CLOSED(0.0, 2.0/3.0), "kappa");
	ADD_PARAMETER(m_fiber, "fiber");
END_FECORE_CLASS();

FEContinuousElasticDamage::FEContinuousElasticDamage(FEModel* fem) : FEElasticMaterial(fem)
{
	m_c1 = 0.0;
	m_k = 0.0;

	m_a1 = 0.0;
	m_a2 = 0.0;
	m_kappa = 0.0;

	m_fiber = vec3d(1, 0, 0);
}

FEMaterialPoint* FEContinuousElasticDamage::CreateMaterialPointData()
{
	return new FEContinuousElasticDamage::Data(new FEElasticMaterialPoint);
}

double FEContinuousElasticDamage::StrainEnergyDensity(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	mat3ds C = pt.RightCauchyGreen();
	mat3ds C2 = C.sqr();
	double I1 = C.tr();
	double J = pt.m_J;
	double J23 = pow(J, 2. / 3.);
	double lnJ = log(J);

	double I1_tilde = I1 / J23;

	// matrix strain-energy
	double Psi_iso = m_c1 * (I1_tilde - 3.0);
	double Psi_vol = 0.5*m_k*(lnJ*lnJ);
	double Psi_matrix = Psi_iso + Psi_vol;

	// fiber strain energy
	vec3d a0 = m_fiber(mp);
	double I4 = a0 * (C*a0);
	double I5 = a0 * (C2*a0);
	double K3 = I1 * I4 - I5;

	double Psi_fiber = m_a1 * pow((m_kappa*I1 + (1.0-3.0*m_kappa/2.0)*K3)-2.0, m_a2);


	// add it all up
	double Psi = Psi_matrix + Psi_fiber;

	return Psi;
}

mat3ds FEContinuousElasticDamage::Stress(FEMaterialPoint& mp)
{
	// calculate the matrix stress
	mat3ds s_matrix = MatrixStress(mp);

	// calculate the fiber stress
	mat3ds s_fiber = FiberStress(mp);

	// add it all up
	mat3ds s = s_matrix + s_fiber;

	// all done
	return s;
}

mat3ds FEContinuousElasticDamage::MatrixStress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	double J = pt.m_J;
	double Jm53 = pow(J, -5. / 3.);

	mat3ds b = pt.LeftCauchyGreen();
	double I1 = b.tr();

	double p = m_k * log(J) / J;

	mat3dd I(1.0);
	mat3ds s_iso = (b - I*(I1/3.0))*(2.0*m_c1 * Jm53);
	mat3ds s_vol = I * p;

	mat3ds s = s_iso + s_vol;
	
	return s;
}

mat3ds FEContinuousElasticDamage::FiberStress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	mat3d F = pt.m_F;
	mat3ds C = pt.RightCauchyGreen();
	mat3ds C2 = C.sqr();
	mat3ds b = pt.LeftCauchyGreen();
	double I1 = C.tr();
	double J = pt.m_J;

	vec3d a0 = m_fiber(mp);

	vec3d a = F * a0;
	a.unit();

	double I4 = a0 * (C *a0);
	double I5 = a0 * (C2*a0);
	double K3 = I1 * I4 - I5;

	double g1 = 1.0 - 3.0*m_kappa / 2.0;
	double g2 = m_kappa * I1 + g1 * K3 - 2.0;
	if (g2 < 0.0) g2 = 0.0;

	double g = 2.0*m_a1*m_a2*pow(g2, m_a2 - 1.0);

	mat3dd I(1.0);
	mat3ds m = dyad(a);
	mat3ds aob = dyads(a, b*a);
	mat3ds t = b * I4 + m * (I1*I4) - aob * I4;
	mat3ds s = (b*m_kappa + t*g1)*(g / J);

	return s;
}

tens4ds FEContinuousElasticDamage::Tangent(FEMaterialPoint& mp)
{
	// calculate matrix tangent
	tens4ds c_matrix = MatrixTangent(mp);

	// calculate fiber tangent
	tens4ds c_fiber = FiberTangent(mp);

	// add it all up
	tens4ds c = c_matrix + c_fiber;

	// all done
	return c;
}

//-----------------------------------------------------------------------------
//! calculate spatial tangent stiffness at material point, using secant method
mat3ds FEContinuousElasticDamage::SecantStress(FEMaterialPoint& mp)
{
	// extract the deformation gradient
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	mat3d F = pt.m_F;
	double J = pt.m_J;
	mat3ds E = pt.Strain();
	mat3dd I(1);
	mat3d FiT = F.transinv();

	// calculate the 2nd P-K stress at the current deformation gradient
	double W = StrainEnergyDensity(mp);

	// create deformation gradient increment
	double eps = 1e-7;
	vec3d e[3];
	e[0] = vec3d(1, 0, 0); e[1] = vec3d(0, 1, 0); e[2] = vec3d(0, 0, 1);
	mat3ds S(0.0);
	for (int k = 0; k < 3; ++k) {
		for (int l = k; l < 3; ++l) {
			// evaluate incremental stress
			mat3d dF = FiT * ((e[k] & e[l]))*(eps*0.5);
			mat3d F1 = F + dF;
			pt.m_F = F1;
			pt.m_J = pt.m_F.det();

			double dW = StrainEnergyDensity(mp) - W;

			// evaluate the secant modulus
			S(k, l) = 2.0 * dW / eps;
		}
	}

	// push from material to spatial frame
	mat3ds s = pt.push_forward(S);

	// restore values
	pt.m_F = F;
	pt.m_J = J;

	// return secant stress
	return s;
}

tens4ds FEContinuousElasticDamage::FiberTangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	mat3d F = pt.m_F;
	mat3ds C = pt.RightCauchyGreen();
	mat3ds C2 = C.sqr();
	mat3ds b = pt.LeftCauchyGreen();
	double I1 = C.tr();
	double J = pt.m_J;

	vec3d a0 = m_fiber(mp);
	vec3d a = F * a0;
	a.unit();

	double I4 = a0 * (C *a0);
	double I5 = a0 * (C2*a0);
	double K3 = I1 * I4 - I5;

	mat3ds m = dyad(a);
	mat3ds aBa = dyads(a, b*a);
	
	mat3dd I(1.0);
	mat3ds t = b * I4 + m * (I1*I4) - aBa * I4;
	tens4ds p = dyad4s(a, b, a)*(I4*2.0);
	tens4ds bom = dyad1s(b, m);
	tens4ds bot = dyad1s(b, t);
	tens4ds txt = dyad1s(t);
	tens4ds bxb = dyad1s(b);

	double g1 = 1.0 - 3.0*m_kappa / 2.0;
	double g2 = m_kappa * I1 + g1 * K3 - 2.0;
	if (g2 < 0.0) g2 = 0.0;

	tens4ds c = (bxb*(m_kappa*m_kappa) + bot*(m_kappa*g1) + txt*(g1*g1))*((m_a2-1.0)*pow(g2, m_a2-2));
	c += (bom*I4 - p)*(g1*pow(g2, m_a2 - 1.0));
	c *= 4.0*m_a1*m_a2 / J;

	return c;
}

tens4ds FEContinuousElasticDamage::MatrixTangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	double J = pt.m_J;
	double Jm53 = pow(J, -5.0 / 3.0);

	mat3ds b = pt.LeftCauchyGreen();
	double I1 = b.tr();

	mat3ds I(mat3dd(1.0));

	double p = m_k * log(J) / J;
	double dp = m_k * (1.0 - log(J))/(J*J);

	tens4ds I4 = dyad4s(I);
	tens4ds IxI = dyad1s(I);

	tens4ds c_iso = (I4*I1 - dyad1s(I, b) + IxI*(I1 / 3.0))*(4.0*m_c1*Jm53 / 3.0);

	tens4ds c_vol = IxI * (dp*J) + (IxI - I4*2.0)*p;

	tens4ds c = c_iso + c_vol;

	return c;
}
