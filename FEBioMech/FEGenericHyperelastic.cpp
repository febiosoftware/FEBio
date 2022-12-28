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
#include "FEGenericHyperelastic.h"
#include <FECore/MMath.h>
#include <FECore/MObj2String.h>
#include <FECore/log.h>

BEGIN_FECORE_CLASS(FEGenericHyperelastic, FEElasticMaterial)
	ADD_PARAMETER(m_exp, "W")->setUnits(UNIT_ENERGY);
END_FECORE_CLASS();

FEGenericHyperelastic::FEGenericHyperelastic(FEModel* fem) : FEElasticMaterial(fem)
{
}

bool FEGenericHyperelastic::Init()
{
	if (BuildMathExpressions() == false) return false;
	return FEElasticMaterial::Init();
}

bool FEGenericHyperelastic::BuildMathExpressions()
{
	vector<string> vars = { "I1", "I2", "J", "X", "Y", "Z"};

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
    MITEM WJ = MSimplify(MDerive(m_W.GetExpression(), *m_W.Variable(2), 1));
	m_W1.AddVariables(vars); m_W1.SetExpression(W1);
	m_W2.AddVariables(vars); m_W2.SetExpression(W2);
	m_WJ.AddVariables(vars); m_WJ.SetExpression(WJ);

    MITEM W11 = MDerive(m_W1.GetExpression(), *m_W1.Variable(0), 1);
    MITEM W12 = MDerive(m_W1.GetExpression(), *m_W1.Variable(1), 1);
    MITEM W22 = MDerive(m_W2.GetExpression(), *m_W2.Variable(1), 1);
	m_W11.AddVariables(vars); m_W11.SetExpression(W11);
	m_W12.AddVariables(vars); m_W12.SetExpression(W12);
	m_W22.AddVariables(vars); m_W22.SetExpression(W22);

    MITEM WJJ = MDerive(m_WJ.GetExpression(), *m_WJ.Variable(2), 1);
	m_WJJ.AddVariables(vars); m_WJJ.SetExpression(WJJ);

#ifdef _DEBUG
	MObj2String o2s;
	string sW1 = o2s.Convert(m_W1); feLog("W1  = %s\n", sW1.c_str());
	string sW2 = o2s.Convert(m_W2); feLog("W2  = %s\n", sW2.c_str());
	string sWJ = o2s.Convert(m_WJ); feLog("WJ  = %s\n", sWJ.c_str());
	string sW11 = o2s.Convert(m_W11); feLog("W11 = %s\n", sW11.c_str());
	string sW12 = o2s.Convert(m_W12); feLog("W12 = %s\n", sW12.c_str());
	string sW22 = o2s.Convert(m_W22); feLog("W22 = %s\n", sW22.c_str());
	string sWJJ = o2s.Convert(m_WJJ); feLog("WJJ = %s\n", sWJJ.c_str());
#endif

	return true;
}


// serialization
void FEGenericHyperelastic::Serialize(DumpStream& ar)
{
	FEElasticMaterial::Serialize(ar);
	if ((ar.IsShallow() == false) && ar.IsLoading())
	{
		bool b = BuildMathExpressions();
		assert(b);
	}
}

mat3ds FEGenericHyperelastic::Stress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	mat3d& F = pt.m_F;
	double J = pt.m_J;

	mat3ds B = pt.LeftCauchyGreen();
	mat3ds B2 = B.sqr();

	double I1 = B.tr();
	double I2 = 0.5*(I1*I1 - B2.tr());

	vector<double> v = { I1, I2, J, mp.m_r0.x, mp.m_r0.y, mp.m_r0.z };
	for (int i = 0; i < m_param.size(); ++i) v.push_back(*m_param[i]);

	double W1 = m_W1.value_s(v);
	double W2 = m_W2.value_s(v);
	double WJ = m_WJ.value_s(v);

	mat3dd I(1.0);

	mat3ds s = (B*(W1 + W2*I1) - B2*W2)*(2.0/J) + WJ*I;

	return s;
}

tens4ds FEGenericHyperelastic::Tangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	mat3d& F = pt.m_F;
	double J = pt.m_J;

	mat3ds B = pt.LeftCauchyGreen();
	mat3ds B2 = B.sqr();

	double I1 = B.tr();
	double I2 = 0.5*(I1*I1 - B2.tr());

	vector<double> v = { I1, I2, J, mp.m_r0.x, mp.m_r0.y, mp.m_r0.z };
	for (int i = 0; i < m_param.size(); ++i) v.push_back(*m_param[i]);

	double W1 = m_W1.value_s(v);
	double W2 = m_W2.value_s(v);
	double WJ = m_WJ.value_s(v);

	double W11 = m_W11.value_s(v);
	double W22 = m_W22.value_s(v);
	double W12 = m_W12.value_s(v);

	double WJJ = m_WJJ.value_s(v);

	mat3dd I(1.0);
	tens4ds IxI = dyad1s(I);

	tens4ds BxB = dyad1s(B);
	tens4ds B2xB2 = dyad1s(B2);

	tens4ds BoB2 = dyad1s(B, B2);

	tens4ds Ib = dyad4s(B);
	tens4ds I4 = dyad4s(I);

	tens4ds cw = BxB*(W11 + 2*W12*I1 + W22*I1*I1 + W2) - BoB2*(W12 + W22*I1) + B2xB2*W22 - Ib*W2;

	tens4ds cp = IxI*(WJJ*J + WJ) - I4*(2 * WJ);

	tens4ds c = cp + cw*(4.0 / J);

	return c;
}

double FEGenericHyperelastic::StrainEnergyDensity(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	mat3d& F = pt.m_F;
	double J = pt.m_J;

	mat3ds B = pt.LeftCauchyGreen();
	mat3ds B2 = B.sqr();

	double I1 = B.tr();
	double I2 = 0.5*(I1*I1 - B2.tr());

	vector<double> v = { I1, I2, J, mp.m_r0.x, mp.m_r0.y, mp.m_r0.z };
	for (int i = 0; i < m_param.size(); ++i) v.push_back(*m_param[i]);

	double W = m_W.value_s(v);

	return W;
}
