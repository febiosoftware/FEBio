/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, The Trustees of Columbia University in 
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
#include "FEFunction1D.h"
#include "DumpStream.h"
#include "MMath.h"
#include "MObj2String.h"

REGISTER_SUPER_CLASS(FEFunction1D, FEFUNCTION1D_ID);

FEFunction1D::FEFunction1D(FEModel* fem) : FECoreBase(fem)
{
}

double FEFunction1D::derive(double x) const
{
	const double eps = 1e-6;
	return (value(x + eps) - value(x))/eps;
}

void FEFunction1D::Serialize(DumpStream& ar)
{
	FECoreBase::Serialize(ar);

	if (ar.IsSaving())
	{
	}
	else
	{
	}
}


//=============================================================================
BEGIN_FECORE_CLASS(FELinearFunction, FEFunction1D)
	ADD_PARAMETER(m_slope, "slope");
	ADD_PARAMETER(m_intercept, "intercept");
END_FECORE_CLASS();

//=============================================================================
BEGIN_FECORE_CLASS(FEMathFunction, FEFunction1D)
	ADD_PARAMETER(m_s, "math");
END_FECORE_CLASS();

FEMathFunction::FEMathFunction(FEModel* fem) : FEFunction1D(fem)
{
	m_s = "0";
}

bool FEMathFunction::Init()
{
	if (m_exp.Create(m_s, "true") == false) return false;
	if (m_exp.Variables() > 1) return false;
	if (m_exp.Variables() < 1) m_exp.AddVariable("x");

	m_dexp.AddVariable(m_exp.Variable(0)->Name());
    if (m_exp.Variables() == 1) {
        MITEM mi = MDerive(m_exp.GetExpression(), *m_exp.Variable(0));
		m_dexp.SetExpression(mi);
    }
	else
		m_dexp.Create("0");

#ifdef _DEBUG
	MObj2String o2s;
	string s = o2s.Convert(m_dexp);
#endif

	return FEFunction1D::Init();
}

FEFunction1D* FEMathFunction::copy()
{
	FEMathFunction* m = new FEMathFunction(GetFEModel());
	m->m_s = m_s;
	m->m_exp = m_exp;
	m->m_dexp = m_dexp;
	return m;
}

double FEMathFunction::value(double t) const
{
	vector<double> var(1, t);
	return m_exp.value_s(var);
}

double FEMathFunction::derive(double t) const
{
	vector<double> var(1, t);
	return m_dexp.value_s(var);
}
