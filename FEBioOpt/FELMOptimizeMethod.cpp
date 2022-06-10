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
#include "FELMOptimizeMethod.h"
#include "FEOptimizeData.h"
#include "FEOptimizeInput.h"
#include "FECore/FEAnalysis.h"
#include "FECore/log.h"

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FELMOptimizeMethod, FEOptimizeMethod)
	ADD_PARAMETER(m_objtol, "obj_tol"     );
	ADD_PARAMETER(m_fdiff , "f_diff_scale");
	ADD_PARAMETER(m_nmax  , "max_iter"    );
	ADD_PARAMETER(m_bcov  , "print_cov"   );
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FELMOptimizeMethod* FELMOptimizeMethod::m_pThis = 0;

// forward declarations
void mrqmin(vector<double>& x, 
			vector<double>& y, 
			vector<double>& sig, 
			vector<double>& a, 
			matrix& covar, 
			matrix& alpha, 
			vector<double>& oneda,
			vector<double>& atry,
			vector<double>& beta,
			vector<double>& da,
			double& chisq,
			void funcs(vector<double>& , vector<double>&, vector<double>&, matrix&),
			double& alamda);

void mrqcof(vector<double>& x, 
			vector<double>& y, 
			vector<double>& sig, 
			vector<double>& a, 
			matrix& alpha,
			vector<double>& beta,
			double& chisq,
			void funcs(vector<double>& , vector<double>&, vector<double>&, matrix&));

//-----------------------------------------------------------------------------
FELMOptimizeMethod::FELMOptimizeMethod(FEModel* fem) : FEOptimizeMethod(fem)
{
	m_objtol = 0.001;
	m_fdiff  = 0.001;
	m_nmax   = 100;
	m_bcov   = 0;
	m_loglevel = LogLevel::LOG_NEVER;
}

//-----------------------------------------------------------------------------
bool FELMOptimizeMethod::Solve(FEOptimizeData *pOpt, vector<double>& amin, vector<double>& ymin, double* minObj)
{
	m_pOpt = pOpt;
	FEOptimizeData& opt = *pOpt;

	// set the variables
	int ma = opt.InputParameters();
	vector<double> a(ma);
	for (int i=0; i<ma; ++i)
	{
		FEInputParameter& var = *opt.GetInputParameter(i);
		a[i] = var.GetValue();
	}

	// set the data
	FEObjectiveFunction& obj = opt.GetObjective();
	int ndata = obj.Measurements();
	vector<double> x(ndata), y(ndata);
	obj.GetMeasurements(y);

	// we don't really need x so we create a dummy array
	for (int i = 0; i<ndata; ++i) x[i] = i;

	// set the sigma's
	// for now we set them all to 1
	vector<double> sig(ndata);
	for (int i = 0; i<ndata; ++i) sig[i] = 1;

	// allocate matrices
	matrix covar(ma, ma), alpha(ma, ma);

	// allocate vectors for mrqmin
	vector<double> oneda(ma), atry(ma), beta(ma), da(ma);

	// set the this pointer
	m_pThis = this;

	opt.m_niter = 0;

	// return value
	double fret = 0.0;

	int niter = 1;

	FEModel* fem = pOpt->GetFEModel();

	try
	{
		// do the first call with lamda to intialize the minimization
		double alamda = -1.0;
		feLogEx(fem, "\n----- Major Iteration: %d -----\n", 0);
		mrqmin(x, y, sig, a, covar, alpha, oneda, atry, beta, da, fret, objfun, alamda);

		// repeat until converged
		double fprev = fret, lam1 = alamda;
		bool bconv = false;
		do
		{
			feLogEx(fem, "\n----- Major Iteration: %d -----\n", niter);
			mrqmin(x, y, sig, a, covar, alpha, oneda, atry, beta, da, fret, objfun, alamda);

			if (alamda < lam1)
			{
				if (niter != 1)
				{
					double df = (fprev - fret)/(fprev + fret + 1);
					if ( df < m_objtol) bconv = true;
					feLogEx(fem, "objective value: %lg (diff = %lg)\n\n", fret, df);
				}
			}
			else feLogEx(fem, "\n objective value: %lg\n\n", fret);

			fprev = fret;
			lam1 = alamda;

			++niter;
		}
		while ((bconv == false) && (niter < m_nmax));

		// do final call with lamda = 0
		alamda = 0.0;
		mrqmin(x, y, sig, a, covar, alpha, oneda, atry, beta, da, fret, objfun, alamda);

	}
	catch (FEErrorTermination)
	{
		feLogErrorEx(fem, "FEBio error terminated. Parameter optimization cannot continue.");
		return false;
	}

	// store optimal values
	amin = a;
	ymin = m_yopt;
	if (minObj) *minObj = fret;

	return true;
}

//-----------------------------------------------------------------------------
void FELMOptimizeMethod::ObjFun(vector<double>& x, vector<double>& a, vector<double>& y, matrix& dyda)
{
	// get the optimization data
	FEOptimizeData& opt = *m_pOpt;

	// poor man's box constraints
	double dir = 1;	// forward difference by default
	for (int i=0; i<opt.InputParameters(); ++i)
	{
		FEInputParameter& var = *opt.GetInputParameter(i);
		if (a[i] < var.MinValue()) {
			a[i] = var.MinValue();
		} else if (a[i] >= var.MaxValue()) {
			a[i] = var.MaxValue();
			dir = -1;	// use backward difference
		}
	}
	
	// evaluate at a
	if (opt.FESolve(a) == false) throw FEErrorTermination();
	
	opt.GetObjective().Evaluate(y);
	m_yopt = y;

	// now calculate the derivatives using forward differences
	int ndata = (int)x.size();
	vector<double> a1(a);
	vector<double> y1(ndata);
	int ma = (int)a.size();
	for (int i=0; i<ma; ++i)
	{
		FEInputParameter& var = *opt.GetInputParameter(i);

		double b = var.ScaleFactor();

		a1[i] = a1[i] + dir*m_fdiff*(fabs(b) + fabs(a[i]));
		assert(a1[i] != a[i]);

		if (opt.FESolve(a1) == false) throw FEErrorTermination();
		opt.GetObjective().Evaluate(y1);
		for (int j=0; j<ndata; ++j) dyda[j][i] = (y1[j] - y[j])/(a1[i] - a[i]);
		a1[i] = a[i];
	}
}

//-----------------------------------------------------------------------------
void mrqmin(vector<double>& x, 
			vector<double>& y, 
			vector<double>& sig, 
			vector<double>& a, 
			matrix& covar, 
			matrix& alpha, 
			vector<double>& oneda,
			vector<double>& atry,
			vector<double>& beta,
			vector<double>& da,
			double& chisq,
			void funcs(vector<double>& , vector<double>&, vector<double>&, matrix&),
			double& alamda)
{
	static double ochisq;
	int j, k, l;

	int ma = (int)a.size();
	if (alamda < 0)
	{
		alamda = 0.001;
		mrqcof(x, y, sig, a, alpha, beta, chisq, funcs);
		ochisq = chisq;
		for (j=0; j<ma; j++) atry[j] = a[j];
	}
	matrix temp(ma, ma);
	for (j=0; j<ma; j++)
	{
		for (k=0; k<ma; k++) covar[j][k] = alpha[j][k];
		covar[j][j] = alpha[j][j]*(1.0 + alamda);
		for (k=0; k<ma; k++) temp[j][k] = covar[j][k];
		oneda[j] = beta[j];
	}
	matrix tempi = temp.inverse();
	oneda = tempi*oneda;
	for (j=0; j<ma; j++)
	{
		for (k=0; k<ma; k++) covar[j][k] = tempi[j][k];
		da[j] = oneda[j];
	}

	if (alamda == 0.0) return;

	for (j=0, l=0; l<ma; l++) atry[l] = a[l] + da[j++];
	mrqcof(x, y, sig, atry, covar, da, chisq, funcs);
	if (chisq < ochisq)
	{
		alamda *= 0.1;
		ochisq = chisq;
		for (j=0; j<ma; j++)
		{
			for (k=0; k<ma; k++) alpha[j][k] = covar[j][k];
			beta[j] = da[j];
		}
		for (l=0; l<ma; l++) a[l] = atry[l];
	}
	else
	{
		alamda *= 10.0;
		chisq = ochisq;
	}
}

//-----------------------------------------------------------------------------
void mrqcof(vector<double>& x, 
			vector<double>& y, 
			vector<double>& sig, 
			vector<double>& a, 
			matrix& alpha,
			vector<double>& beta,
			double& chisq,
			void funcs(vector<double>& , vector<double>&, vector<double>&, matrix&))
{
	int i, j, k, l, m;
	double wt, sig2i, dy;

	int ndata = (int)x.size();
	int ma = (int)a.size();
	for (j=0; j<ma; j++)
	{
		for (k=0; k<=j; k++) alpha[j][k] = 0.0;
		beta[j] = 0.0;
	}

	vector<double> ymod(ndata);
	matrix dyda(ndata, ma);
	funcs(x, a, ymod, dyda);
	
	chisq = 0.0;
	for (i=0; i<ndata; i++)
	{
		sig2i = 1.0 / (sig[i]*sig[i]);
		dy = y[i] - ymod[i];
		for (j=0, l=0; l<ma; l++)
		{
			wt = dyda[i][l]*sig2i;
			for (k=0, m=0; m<l+1; m++) alpha[j][k++] += wt*dyda[i][m];
			beta[j++] += dy*wt;
		}
		chisq += dy*dy*sig2i;
	}
	for (j=1; j<ma; j++)
		for (k=0; k<j; k++) alpha[k][j] = alpha[j][k];
}
