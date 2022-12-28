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
#include "FEGenericTransIsoHyperelasticUC.h"
#include <FECore/MMath.h>
#include <FECore/MObj2String.h>
#include <FECore/log.h>
#include <FECore/FEConstValueVec3.h>

BEGIN_FECORE_CLASS(FEGenericTransIsoHyperelasticUC, FEUncoupledMaterial)
	ADD_PARAMETER(m_exp, "W");
	ADD_PARAMETER(m_printDerivs, "print_derivs");

	ADD_PROPERTY(m_fiber, "fiber");
END_FECORE_CLASS();

FEGenericTransIsoHyperelasticUC::FEGenericTransIsoHyperelasticUC(FEModel* fem) : FEUncoupledMaterial(fem)
{
	m_fiber = nullptr;
	m_printDerivs = false;
}

bool FEGenericTransIsoHyperelasticUC::Init()
{
	vector<string> vars = { "I1", "I2", "I4", "I5" };

	// add all user parameters
	FEParameterList& pl = GetParameterList();
	FEParamIterator pi = pl.first();
	m_param.clear();
	for (int i = 0; i < pl.Parameters(); ++i, ++pi)
	{
		FEParam& p = *pi;
		if (p.GetFlags() & FEParamFlag::FE_PARAM_USER)
		{
			vars.push_back(p.name());
			m_param.push_back((double*)p.data_ptr());
		}
	}

	// create math object
	m_W.AddVariables(vars);
	if (m_W.Create(m_exp) == false) return false;

	// calculate all derivatives
	MITEM W1 = MSimplify(MDerive(m_W.GetExpression(), *m_W.Variable(0), 1));
	MITEM W2 = MSimplify(MDerive(m_W.GetExpression(), *m_W.Variable(1), 1));
	MITEM W4 = MSimplify(MDerive(m_W.GetExpression(), *m_W.Variable(2), 1));
	MITEM W5 = MSimplify(MDerive(m_W.GetExpression(), *m_W.Variable(3), 1));
	MITEM WJ = MSimplify(MDerive(m_W.GetExpression(), *m_W.Variable(4), 1));
	m_W1.AddVariables(vars); m_W1.SetExpression(W1);
	m_W2.AddVariables(vars); m_W2.SetExpression(W2);
	m_W4.AddVariables(vars); m_W4.SetExpression(W4);
	m_W5.AddVariables(vars); m_W5.SetExpression(W5);

	MITEM W11 = MDerive(m_W1.GetExpression(), *m_W1.Variable(0), 1);
	MITEM W12 = MDerive(m_W1.GetExpression(), *m_W1.Variable(1), 1);
	MITEM W14 = MDerive(m_W1.GetExpression(), *m_W1.Variable(2), 1);
	MITEM W15 = MDerive(m_W1.GetExpression(), *m_W1.Variable(3), 1);
	MITEM W22 = MDerive(m_W2.GetExpression(), *m_W2.Variable(1), 1);
	MITEM W24 = MDerive(m_W2.GetExpression(), *m_W2.Variable(2), 1);
	MITEM W25 = MDerive(m_W2.GetExpression(), *m_W2.Variable(3), 1);
	MITEM W44 = MDerive(m_W4.GetExpression(), *m_W4.Variable(2), 1);
	MITEM W45 = MDerive(m_W4.GetExpression(), *m_W4.Variable(3), 1);
	MITEM W55 = MDerive(m_W5.GetExpression(), *m_W5.Variable(3), 1);
	m_W11.AddVariables(vars); m_W11.SetExpression(W11);
	m_W12.AddVariables(vars); m_W12.SetExpression(W12);
	m_W14.AddVariables(vars); m_W14.SetExpression(W14);
	m_W15.AddVariables(vars); m_W15.SetExpression(W15);
	m_W22.AddVariables(vars); m_W22.SetExpression(W22);
	m_W24.AddVariables(vars); m_W24.SetExpression(W24);
	m_W25.AddVariables(vars); m_W25.SetExpression(W25);
	m_W44.AddVariables(vars); m_W44.SetExpression(W44);
	m_W45.AddVariables(vars); m_W45.SetExpression(W45);
	m_W55.AddVariables(vars); m_W55.SetExpression(W55);

	if (m_printDerivs)
	{
		feLog("\nStrain energy and derivatives for material %d (%s):\n", GetID(), GetName().c_str());
		MObj2String o2s;
		string sW = o2s.Convert(m_W); feLog("W = %s\n", sW.c_str());
		string sW1 = o2s.Convert(m_W1); feLog("W1  = %s\n", sW1.c_str());
		string sW2 = o2s.Convert(m_W2); feLog("W2  = %s\n", sW2.c_str());
		string sW4 = o2s.Convert(m_W4); feLog("W4  = %s\n", sW4.c_str());
		string sW5 = o2s.Convert(m_W5); feLog("W5  = %s\n", sW5.c_str());
		string sW11 = o2s.Convert(m_W11); feLog("W11 = %s\n", sW11.c_str());
		string sW12 = o2s.Convert(m_W12); feLog("W12 = %s\n", sW12.c_str());
		string sW14 = o2s.Convert(m_W14); feLog("W14 = %s\n", sW14.c_str());
		string sW15 = o2s.Convert(m_W15); feLog("W15 = %s\n", sW15.c_str());
		string sW22 = o2s.Convert(m_W22); feLog("W22 = %s\n", sW22.c_str());
		string sW24 = o2s.Convert(m_W24); feLog("W24 = %s\n", sW24.c_str());
		string sW25 = o2s.Convert(m_W25); feLog("W25 = %s\n", sW25.c_str());
		string sW44 = o2s.Convert(m_W44); feLog("W44 = %s\n", sW44.c_str());
		string sW45 = o2s.Convert(m_W45); feLog("W45 = %s\n", sW45.c_str());
		string sW55 = o2s.Convert(m_W55); feLog("W55 = %s\n", sW55.c_str());
	}

	return FEElasticMaterial::Init();
}

mat3ds FEGenericTransIsoHyperelasticUC::DevStress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	mat3d& F = pt.m_F;
	double J = pt.m_J;

	mat3ds B = pt.DevLeftCauchyGreen();
	mat3ds B2 = B.sqr();

	// get the material fiber axis
	vec3d a0 = m_fiber->unitVector(mp);

	// get the spatial fiber axis
	vec3d a = pt.m_F*a0;
	double lam = a.unit()*pow(J, -1./3.);

	// evaluate the invariants
	double I1 = B.tr();
	double I2 = 0.5*(I1*I1 - B2.tr());
	double I4 = lam*lam;
	double I5 = I4*(a*(B*a));

	// create the parameter list
	vector<double> v = { I1, I2, I4, I5 };
	for (int i = 0; i < m_param.size(); ++i) v.push_back(*m_param[i]);

	// evaluate the strain energy derivatives
	double W1 = m_W1.value_s(v);
	double W2 = m_W2.value_s(v);
	double W4 = m_W4.value_s(v);
	double W5 = m_W5.value_s(v);

	mat3dd I(1.0);
	mat3ds AxA = dyad(a);
	mat3ds aBa = dyads(a, B*a);

	mat3ds T = (B*(W1 + W2*I1) - B2*W2 + AxA*(W4*I4) + aBa*(W5*I4));

	mat3ds s = T.dev()*(2.0 / J);

	// all done!
	return s;
}

tens4ds FEGenericTransIsoHyperelasticUC::DevTangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	mat3d& F = pt.m_F;
	double J = pt.m_J;

	mat3ds B = pt.DevLeftCauchyGreen();
	mat3ds B2 = B.sqr();

	// get the material fiber axis
	vec3d a0 = m_fiber->unitVector(mp);

	// get the spatial fiber axis
	vec3d a = pt.m_F*a0;
	double lam = a.unit()*pow(J, -1./3.);

	// evaluate the invariants
	double I1 = B.tr();
	double I2 = 0.5*(I1*I1 - B2.tr());
	double I4 = lam*lam;
	double I5 = I4*(a*(B*a));

	// evaluate parameters
	vector<double> v = { I1, I2, I4, I5 };
	for (int i = 0; i < m_param.size(); ++i) v.push_back(*m_param[i]);

	// evaluate strain energy derivatives
	double W1 = m_W1.value_s(v);
	double W2 = m_W2.value_s(v);
	double W4 = m_W4.value_s(v);
	double W5 = m_W5.value_s(v);

	double W11 = m_W11.value_s(v);
	double W12 = m_W12.value_s(v);
	double W14 = m_W14.value_s(v);
	double W15 = m_W15.value_s(v);
	double W22 = m_W22.value_s(v);
	double W24 = m_W24.value_s(v);
	double W25 = m_W25.value_s(v);
	double W44 = m_W44.value_s(v);
	double W45 = m_W45.value_s(v);
	double W55 = m_W55.value_s(v);

	// a few tensors we'll need
	mat3dd I(1.0);
	tens4ds IxI = dyad1s(I);
	tens4ds BxB = dyad1s(B);
	tens4ds B2xB2 = dyad1s(B2);
	tens4ds BoB2 = dyad1s(B, B2);
	tens4ds Ib = dyad4s(B);
	tens4ds II = dyad4s(I);
	mat3ds A = dyad(a);
	mat3ds aBa = dyads(a, B*a);
	tens4ds Baa = dyad1s(B, A);
	tens4ds BaBa = dyad1s(B, aBa);
	tens4ds B2A = dyad1s(B2, A);
	tens4ds B2aBa = dyad1s(B2, aBa);
	tens4ds AxA = dyad1s(A);
	tens4ds AaBa = dyad1s(A, aBa);
	tens4ds aBaaBa = dyad1s(aBa, aBa);
	tens4ds AB = dyad5s(A, B);

	// evaluate deviatoric stress
	mat3ds T = (B*(W1 + W2*I1) - B2*W2 + A*(W4*I4) + aBa*(W5*I4))*(2./J);
	mat3ds devT = T.dev();

	// stiffness contribution from isotropic terms
	tens4ds dw_iso = BxB*(W11 + 2.0*W12*I1 + W22*I1*I1 + W2) - BoB2*(W12 + W22*I1) + B2xB2*W22 - Ib*W2;

	// stiffness constribution from anisotropic terms 
	tens4ds dw_ani = Baa*(I4*(W14 + W24*I1)) \
		+ BaBa*(I4*(W15 + W25*I1)) \
		- B2A*(I4*W24) \
		- B2aBa*(I4*W25) \
		+ AxA*(I4*I4*W44) \
		+ AaBa*(I4*I4*W45) \
		+ aBaaBa*(I4*I4*W55) \
		+ AB*(I4*W5);

	// isochoric contribution
	tens4ds dw = dw_iso + dw_ani;
	mat3ds dwoI = dw.contract();
	double trDw = dw.tr();
	tens4ds cw = dw - dyad1s(dwoI, I) / 3.0 + IxI*(trDw / 9.0);

	// put it all together
	tens4ds c = (II - IxI/3.0)*(2.0/3.0*T.tr()) - dyad1s(devT, I)*(2.0/3.0) + cw*(4. / J);

	// all done
	return c;
}

double FEGenericTransIsoHyperelasticUC::DevStrainEnergyDensity(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	mat3d& F = pt.m_F;
	double J = pt.m_J;

	mat3ds B = pt.DevLeftCauchyGreen();
	mat3ds B2 = B.sqr();

	// get the material fiber axis
	vec3d a0 = m_fiber->unitVector(mp);

	// get the spatial fiber axis
	vec3d a = pt.m_F*a0;
	double lam = a.unit()*pow(J, -1./3.);

	// evaluate the invariants
	double I1 = B.tr();
	double I2 = 0.5*(I1*I1 - B2.tr());
	double I4 = lam*lam;
	double I5 = I4*(a*(B*a));

	// evaluate parameters
	vector<double> v = { I1, I2, I4, I5 };
	for (int i = 0; i < m_param.size(); ++i) v.push_back(*m_param[i]);

	double W = m_W.value_s(v);

	return W;
}
