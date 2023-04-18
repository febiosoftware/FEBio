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
#include "FEGenericHyperelasticUC.h"
#include <FECore/MMath.h>
#include <FECore/MObj2String.h>
#include <FECore/log.h>

BEGIN_FECORE_CLASS(FEGenericHyperelasticUC, FEUncoupledMaterial)
	ADD_PARAMETER(m_exp, "W")->setUnits(UNIT_ENERGY);
	ADD_PARAMETER(m_printDerivs, "print_derivs");
END_FECORE_CLASS();

FEGenericHyperelasticUC::FEGenericHyperelasticUC(FEModel* fem) : FEUncoupledMaterial(fem)
{
	m_printDerivs = false;
}

// initialization
bool FEGenericHyperelasticUC::Init()
{
	if (BuildMathExpressions() == false) return false;
	return FEUncoupledMaterial::Init();
}

bool FEGenericHyperelasticUC::BuildMathExpressions()
{
	vector<string> vars = { "I1", "I2", "X", "Y", "Z" };

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

	m_W.AddVariables(vars);
	if (m_W.Create(m_exp) == false) return false;

    MITEM W1 = MSimplify(MDerive(m_W.GetExpression(), *m_W.Variable(0), 1));
    MITEM W2 = MSimplify(MDerive(m_W.GetExpression(), *m_W.Variable(1), 1));
	m_W1.AddVariables(vars); m_W1.SetExpression(W1);
	m_W2.AddVariables(vars); m_W2.SetExpression(W2);

    MITEM W11 = MDerive(m_W1.GetExpression(), *m_W1.Variable(0), 1);
    MITEM W12 = MDerive(m_W1.GetExpression(), *m_W1.Variable(1), 1);
    MITEM W22 = MDerive(m_W2.GetExpression(), *m_W2.Variable(1), 1);
	m_W11.AddVariables(vars); m_W11.SetExpression(W11);
	m_W12.AddVariables(vars); m_W12.SetExpression(W12);
	m_W22.AddVariables(vars); m_W22.SetExpression(W22);

	if (m_printDerivs)
	{
		feLog("\nStrain energy and derivatives for material %d (%s):\n", GetID(), GetName().c_str());
		MObj2String o2s;
		string sW = o2s.Convert(m_W); feLog("W = %s\n", sW.c_str());
		string sW1 = o2s.Convert(m_W1); feLog("W1  = %s\n", sW1.c_str());
		string sW2 = o2s.Convert(m_W2); feLog("W2  = %s\n", sW2.c_str());
		string sW11 = o2s.Convert(m_W11); feLog("W11 = %s\n", sW11.c_str());
		string sW12 = o2s.Convert(m_W12); feLog("W12 = %s\n", sW12.c_str());
		string sW22 = o2s.Convert(m_W22); feLog("W22 = %s\n", sW22.c_str());
	}

	return true;
}

// serialization
void FEGenericHyperelasticUC::Serialize(DumpStream& ar)
{
	FEUncoupledMaterial::Serialize(ar);
	if ((ar.IsShallow() == false) && ar.IsLoading())
	{
		bool b = BuildMathExpressions();
		assert(b);
	}
}

//! Deviatoric Cauchy stress
mat3ds FEGenericHyperelasticUC::DevStress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// determinant of deformation gradient
	double J = pt.m_J;

	// calculate deviatoric left Cauchy-Green tensor
	mat3ds B = pt.DevLeftCauchyGreen();

	// calculate square of B
	mat3ds B2 = B.sqr();

	// Invariants of B (= invariants of C)
	// Note that these are the invariants of Btilde, not of B!
	double I1 = B.tr();
	double I2 = 0.5*(I1*I1 - B2.tr());

	// get strain energy derivatives
	vector<double> v = { I1, I2, mp.m_r0.x, mp.m_r0.y, mp.m_r0.z };
	for (int i = 0; i < m_param.size(); ++i) v.push_back(*m_param[i]);
	double W1 = m_W1.value_s(v);
	double W2 = m_W2.value_s(v);

	// calculate T = F*dW/dC*Ft
	mat3ds T = B*(W1 + W2*I1) - B2*W2;

	// return deviatoric Cauchy stress 
	return T.dev()*(2.0 / J);
}

//! Deviatoric spatial Tangent
tens4ds FEGenericHyperelasticUC::DevTangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// determinant of deformation gradient
	double J = pt.m_J;
	double Ji = 1.0 / J;

	// calculate deviatoric left Cauchy-Green tensor
	mat3ds B = pt.DevLeftCauchyGreen();

	// calculate square of B
	mat3ds B2 = B.sqr();

	// Invariants of B (= invariants of C)
	// Note that these are the invariants of Btilde, not of B!
	double I1 = B.tr();
	double I2 = 0.5*(I1*I1 - B2.tr());

	// get strain energy derivatives
	vector<double> v = { I1, I2, mp.m_r0.x, mp.m_r0.y, mp.m_r0.z };
	for (int i = 0; i < m_param.size(); ++i) v.push_back(*m_param[i]);
	double W1 = m_W1.value_s(v);
	double W2 = m_W2.value_s(v);

	double W11 = m_W11.value_s(v);
	double W22 = m_W22.value_s(v);
	double W12 = m_W12.value_s(v);

	// define fourth-order tensors
	mat3dd I(1.0);
	tens4ds IxI = dyad1s(I);
	tens4ds I4 = dyad4s(I);
	tens4ds BxB = dyad1s(B);
	tens4ds B2xB2 = dyad1s(B2);
	tens4ds BoB2 = dyad1s(B, B2);
	tens4ds Ib = dyad4s(B);

	// deviatoric stress
	mat3ds T = (B*(W1 + W2*I1) - B2*W2)*(2.0/J);
	mat3ds devT = T.dev();

	// isochoric stiffness contribution
	tens4ds dw = BxB*(W11 + 2.0*W12*I1 + W22*I1*I1 + W2) - BoB2*(W12 + W22*I1) + B2xB2*W22 - Ib*W2;
	double trDw = dw.tr();
	mat3ds dwoI = dw.contract();
	tens4ds cw = dw - dyad1s(dwoI, I) / 3.0 + IxI*(trDw / 9.0);

	// put it all together
	tens4ds c = (I4 - IxI/3.0)*(2.0*T.tr()/3.0) - dyad1s(devT, I)*(2.0/3.0) + cw*(4.0/J);

	// all done
	return c;
}

//! Deviatoric strain energy density
double FEGenericHyperelasticUC::DevStrainEnergyDensity(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// determinant of deformation gradient
	double J = pt.m_J;

	// calculate deviatoric left Cauchy-Green tensor
	mat3ds B = pt.DevLeftCauchyGreen();

	// calculate square of B
	mat3ds B2 = B.sqr();

	// Invariants of B (= invariants of C)
	// Note that these are the invariants of Btilde, not of B!
	double I1 = B.tr();
	double I2 = 0.5*(I1*I1 - B2.tr());

	// evaluate (deviatoric) strain energy
	vector<double> v = { I1, I2, mp.m_r0.x, mp.m_r0.y, mp.m_r0.z };
	for (int i = 0; i < m_param.size(); ++i) v.push_back(*m_param[i]);
	double W = m_W.value_s(v);

	return W;
}
