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
#include "FEFunction1D.h"
#include "FEModel.h"
#include "DumpStream.h"
#include "MMath.h"
#include "MObj2String.h"
#include "log.h"

FEFunction1D::FEFunction1D(FEModel* fem) : FECoreBase(fem)
{
}

double FEFunction1D::derive(double x) const
{
	const double eps = 1e-6;
	return (value(x + eps) - value(x))/eps;
}

double FEFunction1D::integrate(double a, double b) const
{
	return (b-a)*((value(a) + value(b))/2);
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

// Invert the function f(x) = f0: given f0, return x.
// On input, x is the initial guess.
// On output, x is the solution.
// Function returns true if solution has converged, false otherwise.
// Uses Newton's method.
bool FEFunction1D::invert(const double f0, double &x)
{
    const int maxiter = 100;
    const double errabs = 1e-15;
    const double errrel = 1e-6;
    
    bool cnvgd = false;
    bool done = false;
    int iter = 0;
    double dx;
    while (!done) {
        double f = value(x) - f0;
        double df = derive(x);
        // check if slope is zero
        if (fabs(df) <= errabs) {
            if (fabs(f) > errabs) {
                // if function not zero, this is a local minimum
                done = true;
            }
            else {
                // if function is zero, this is a valid solution
                done = cnvgd = true;
            }
        }
        // if slope is not zero proceed with iterative scheme
        else {
            dx = -f/df;
            x += dx;
            // check for relative or absolute convergence
            if ((fabs(dx) <= errrel*fabs(x)) || (fabs(f) <= errabs)) {
                done = cnvgd = true;
            }
        }
        ++iter;
        if (iter > maxiter) done = true;
    }
    
    return cnvgd;
}

//=============================================================================
BEGIN_FECORE_CLASS(FEConstFunction, FEFunction1D)
	ADD_PARAMETER(m_value, "value");
END_FECORE_CLASS();

//=============================================================================
BEGIN_FECORE_CLASS(FELinearFunction, FEFunction1D)
	ADD_PARAMETER(m_slope, "slope");
	ADD_PARAMETER(m_intercept, "intercept");
END_FECORE_CLASS();

//=============================================================================
BEGIN_FECORE_CLASS(FEStepFunction, FEFunction1D)
	ADD_PARAMETER(m_x0, "x0");
	ADD_PARAMETER(m_leftVal , "left_val");
	ADD_PARAMETER(m_rightVal, "right_val");
END_FECORE_CLASS();

//=============================================================================
BEGIN_FECORE_CLASS(FEMathFunction, FEFunction1D)
	ADD_PARAMETER(m_s, "math");
END_FECORE_CLASS();

FEMathFunction::FEMathFunction(FEModel* fem) : FEFunction1D(fem)
{
	m_s = "0";
}

void FEMathFunction::SetMathString(const std::string& s)
{
	m_s = s;
}

bool FEMathFunction::Init()
{
	if (BuildMathExpressions() == false) return false;
	return FEFunction1D::Init();
}

bool FEMathFunction::BuildMathExpressions()
{
	// process the string
	m_exp.Clear();
	if (m_exp.Create(m_s, true) == false) return false;

	// match variables to model parameters.
	m_var.clear();
	m_ix = -1;
	FEModel* fem = GetFEModel();
	for (int i=0; i<m_exp.Variables(); ++i)
	{
		MVariable* v = m_exp.Variable(i);
		
		ParamString ps(v->Name().c_str());
		FEParamValue param = fem->GetParameterValue(ps);
		if (param.isValid() == false)
		{
			// let's assume this is the independent parameter
			if (m_ix == -1)
			{
				// push a dummy param value
				m_var.push_back(FEParamValue());
				m_ix = i;
			}
			else return false;
		}
		else
		{
			if (param.type() != FE_PARAM_DOUBLE) return false;
			m_var.push_back(param);
		}
	}

	// copy variables to derived expressions
	m_dexp.Clear();
	m_d2exp.Clear();
	for (int i = 0; i < m_exp.Variables(); ++i)
	{
		m_dexp.AddVariable(m_exp.Variable(i)->Name());
		m_d2exp.AddVariable(m_exp.Variable(i)->Name());
	}

	// evaluate first derivative
    if (m_ix != -1) {
        MITEM mi = MDerive(m_exp.GetExpression(), *m_exp.Variable(m_ix));
		m_dexp.SetExpression(mi);
    }
	else
		m_dexp.Create("0");

	// evaluate second derivative
    if (m_ix != -1) {
        MITEM mi = MDerive(m_dexp.GetExpression(), *m_dexp.Variable(m_ix));
        m_d2exp.SetExpression(mi);
    }
    else
        m_d2exp.Create("0");

	return true;
}

void FEMathFunction::Serialize(DumpStream& ar)
{
	FEFunction1D::Serialize(ar);
	if ((ar.IsShallow() == false) && (ar.IsLoading()))
	{
		bool b = BuildMathExpressions();
		assert(b);
	}
}

FEFunction1D* FEMathFunction::copy()
{
	FEMathFunction* m = new FEMathFunction(GetFEModel());
	m->m_s = m_s;
	m->m_ix = m_ix;
	m->m_var = m_var;

	m->m_exp = m_exp;
	m->m_dexp = m_dexp;
    m->m_d2exp = m_d2exp;
	return m;
}

void FEMathFunction::evalParams(std::vector<double>& val, double t) const
{
	val.resize(m_var.size());
	for (int i = 0; i < m_var.size(); ++i)
	{
		if (i == m_ix) val[i] = t;
		else val[i] = m_var[i].value<double>();
	}
}

double FEMathFunction::value(double t) const
{
	vector<double> v;
	evalParams(v, t);
	return m_exp.value_s(v);
}

double FEMathFunction::derive(double t) const
{
	vector<double> v;
	evalParams(v, t);
	return m_dexp.value_s(v);
}

double FEMathFunction::deriv2(double t) const
{
	vector<double> v;
	evalParams(v, t);
	return m_d2exp.value_s(v);
}

