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
#include "FEConstrainedLMOptimizeMethod.h"
#include "FEOptimizeData.h"
#include "FEOptimizeInput.h"
#include "FECore/FEAnalysis.h"
#include "FECore/log.h"

#ifdef HAVE_LEVMAR
#include "levmar.h"

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEConstrainedLMOptimizeMethod, FEOptimizeMethod)
	ADD_PARAMETER(m_objtol, "obj_tol"     );
	ADD_PARAMETER(m_tau   , "tau"         );
	ADD_PARAMETER(m_fdiff , "f_diff_scale");
	ADD_PARAMETER(m_nmax  , "max_iter"    );
	ADD_PARAMETER(m_scaleParams, "scale_parameters");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEConstrainedLMOptimizeMethod::FEConstrainedLMOptimizeMethod(FEModel* fem) : FEOptimizeMethod(fem)
{
	m_tau = 1e-3;
	m_objtol = 0.001;
	m_fdiff  = 0.001;
	m_nmax   = 100;
	m_scaleParams = false;
    m_loglevel = LogLevel::LOG_NEVER;
}

//-----------------------------------------------------------------------------
bool FEConstrainedLMOptimizeMethod::Solve(FEOptimizeData *pOpt, vector<double>& amin, vector<double>& ymin, double* minObj)
{
	m_pOpt = pOpt;
	FEOptimizeData& opt = *pOpt;

	// get the data
	FEObjectiveFunction& obj = opt.GetObjective();
	int ndata = obj.Measurements();
	vector<double> y(ndata, 0);
	obj.GetMeasurements(y);

	// allocate matrices
	int ma = opt.InputParameters();
	matrix covar(ma, ma), alpha(ma, ma);

	opt.m_niter = 0;

	// return value
	double fret = 0.0;

	int niter = 1;

	// if parameter scaling is not used, just set all scale factors to one
	// to retain backward compatibility
	if (m_scaleParams == false)
	{
		for (int i = 0; i < ma; ++i)
		{
			opt.GetInputParameter(i)->ScaleFactor() = 1.0;
		}
	}

	try
	{
		vector<double> p(ma);
		vector<double> s(ma);
		vector<double> lb(ma);
		vector<double> ub(ma);
		for (int i = 0; i < ma; ++i)
		{
			FEInputParameter& v = *opt.GetInputParameter(i);
			double a = v.GetValue();
			double sf = v.ScaleFactor();

			s[i] = sf;
			p[i] = a / sf;

			lb[i] = v.MinValue() / sf;
			ub[i] = v.MaxValue() / sf;
		}

		vector<double> q(ndata);
		for (int i=0; i<ndata; ++i) q[i] = y[i];

		const double tol = m_objtol;
		double opts[5] = {m_tau, tol, tol, tol, m_fdiff};

		int itmax = m_nmax;
		if (opt.Constraints() > 0)
		{
			int NC = opt.Constraints();
			vector<double> A(NC * ma);
			vector<double> b(NC);
			for (int i=0; i<NC; ++i)
			{
				OPT_LIN_CONSTRAINT& con = opt.Constraint(i);
				for (int j=0; j<ma; ++j) A[i*ma + j] = con.a[j] * s[j];
				b[i] = con.b;
			}

			int ret = dlevmar_blec_dif(objfun, p.data(), q.data(), ma, ndata, lb.data(), ub.data(), A.data(), b.data(), NC, 0, itmax, opts, 0, 0, 0, (void*) this);
		}
		else
		{
			int ret = dlevmar_bc_dif(objfun, p.data(), q.data(), ma, ndata, lb.data(), ub.data(), 0, itmax, opts, 0, 0, 0, (void*) this);
		}

		amin.resize(ma);
		for (int i = 0; i < ma; ++i)
		{
			amin[i] = p[i] * s[i];
		}

		// store the optimal values
		fret = obj.Evaluate(m_yopt);
	}
	catch (FEErrorTermination)
	{
		FEModel* fem = pOpt->GetFEModel();
		feLogErrorEx(fem, "FEBio error terminated. Parameter optimization cannot continue.");
		return false;
	}

	// return optimal values
	ymin = m_yopt;
	if (minObj) *minObj = fret;

	return true;
}

//-----------------------------------------------------------------------------
void FEConstrainedLMOptimizeMethod::ObjFun(double* p, double* hx, int m, int n)
{
	// get the optimization data
	FEOptimizeData& opt = *GetOptimizeData();
	FEObjectiveFunction& obj = opt.GetObjective();

	// evaluate at a
	vector<double> a(m);
	for (int i = 0; i < m; ++i)
	{
		FEInputParameter& var = *opt.GetInputParameter(i);
		a[i] = p[i]*var.ScaleFactor();
	}

	// poor man's box constraints
	for (int i = 0; i < opt.InputParameters(); ++i)
	{
		FEInputParameter& var = *opt.GetInputParameter(i);
		if (a[i] < var.MinValue()) {
			feLogEx(opt.GetFEModel(), "Warning: clamping %s to min (was %lg)\n", var.GetName().c_str(), a[i]);
			a[i] = var.MinValue();
		}
		else if (a[i] >= var.MaxValue()) {
			feLogEx(opt.GetFEModel(), "Warning: clamping %s to max (was %lg)\n", var.GetName().c_str(), a[i]);
			a[i] = var.MaxValue();
		}
	}

	// solve the problem
	if (opt.FESolve(a) == false) throw FEErrorTermination();

	// store the measurement vector
	vector<double> y(n, 0.0);
	opt.GetObjective().Evaluate(y);
	for (int i = 0; i < n; ++i) hx[i] = y[i];

	// store the last calculated values
	m_yopt = y;
}

#endif
