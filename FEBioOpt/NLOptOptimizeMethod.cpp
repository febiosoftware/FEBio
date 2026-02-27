/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2026 University of Utah, The Trustees of Columbia University in
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
#include "NLOptOptimizeMethod.h"
#ifdef HAVE_NLOPT
#include <FECore/FEModel.h>
#include <FECore/FEAnalysis.h>
#include <FECore/FESolver.h>
#include "FEOptimizeData.h"
#include "FEOptimizeInput.h" // for FEErrorTermination
#include <nlopt.h>
#endif

BEGIN_FECORE_CLASS(NLOptOptimizeMethod, FEOptimizeMethod)
	ADD_PARAMETER(m_max_iter, "max_iter");
	ADD_PARAMETER(m_ftol_abs, "ftol_abs");
	ADD_PARAMETER(m_ftol_rel, "ftol_rel");
END_FECORE_CLASS();

#ifdef HAVE_NLOPT
static double myfunc(unsigned int n, const double* x, double* grad, void* my_func_data)
{
	NLOptOptimizeMethod* opt = (NLOptOptimizeMethod*)my_func_data;
	return opt->EvalObjective(n, x, grad);
}

NLOptOptimizeMethod::NLOptOptimizeMethod(FEModel* fem, int nlopt_algor) : FEOptimizeMethod(fem), m_nlopt_algorithm(nlopt_algor) {}

bool NLOptOptimizeMethod::Solve(FEOptimizeData* optData, vector<double>& amin, vector<double>& ymin, double* minObj)
{
	if (m_nlopt_algorithm == -1) return false;

	m_optData = optData;
	int ma = optData->InputParameters();

	nlopt_opt opt;
	opt = nlopt_create((nlopt_algorithm)m_nlopt_algorithm, ma); // algorithm and dimensionality

	// set nlopt options
	vector<double> a(ma), al(ma), au(ma);
	for (int i = 0; i < ma; ++i)
	{
		FEInputParameter* var = optData->GetInputParameter(i);
		a[i] = var->InitValue();
		al[i] = var->MinValue();
		au[i] = var->MaxValue();
	}
	nlopt_set_lower_bounds(opt, al.data());
	nlopt_set_upper_bounds(opt, au.data());

	nlopt_set_min_objective(opt, myfunc, this);
	nlopt_set_maxeval(opt, m_max_iter);
	nlopt_set_param(opt, "verbosity", m_print_level);

	//	nlopt_set_xtol_rel(opt, 1e-9);
	if (m_ftol_abs != 0) nlopt_set_ftol_abs(opt, m_ftol_abs);
	if (m_ftol_rel != 0) nlopt_set_ftol_rel(opt, m_ftol_rel);

	double minf; //the minimum objective value upon return
	int result = 0;
	if ((result = nlopt_optimize(opt, a.data(), &minf)) < 0) {
		printf("nlopt failed!\n");
	}
	else {
		printf("nlopt finished!\n");

		FEObjectiveFunction& obj = optData->GetObjective();
		double f = obj.Evaluate(ymin);
		amin = a;
	}

	if (minObj) *minObj = minf;

	nlopt_destroy(opt);
	return (result >= 0);
}

double NLOptOptimizeMethod::EvalObjective(unsigned int n, const double* x, double* grad)
{
	FEOptimizeData& opt = *m_optData;
	FEObjectiveFunction& obj = opt.GetObjective();

	// evaluate at a
	vector<double> a(n);
	for (int i = 0; i < n; ++i)
	{
		FEInputParameter& var = *opt.GetInputParameter(i);
		a[i] = x[i] * var.ScaleFactor();
	}

	// run the FEBio model with the new parameters
	if (opt.FESolve(a) == false) throw FEErrorTermination();

	// evaluate objective
	double f = obj.Evaluate();

	// calculate a finite difference approximation of the gradient if requested
	if (grad) {
		EvalApproxGradient(n, x, grad);
	}

	return f;
}

void NLOptOptimizeMethod::EvalApproxGradient(unsigned int n, const double* x, double* grad)
{
	FEOptimizeData& opt = *m_optData;
	FEObjectiveFunction& obj = opt.GetObjective();

	const double eps = 1e-6;
	vector<double> x1(n);
	vector<double> x2(n);
	for (int i = 0; i < n; ++i)
	{
		// perturb up
		for (int j = 0; j < n; ++j) x1[j] = x[j];
		x1[i] += eps;
		// run the forward model
		bool b = opt.FESolve(x1);
		double f1 = obj.Evaluate();

		// perturb down
		for (int j = 0; j < n; ++j) x2[j] = x[j];
		x2[i] -= eps;

		// run the forward model
		b = opt.FESolve(x2);
		double f2 = obj.Evaluate();
		double dfdx = (f1 - f2) / (2 * eps);
		grad[i] += dfdx;
	}
}

FEMMAOptimizeMethod::FEMMAOptimizeMethod(FEModel* fem) : NLOptOptimizeMethod(fem, NLOPT_LD_MMA) {}
FENelderMeadOptMethod::FENelderMeadOptMethod(FEModel* fem) : NLOptOptimizeMethod(fem, NLOPT_LN_NELDERMEAD) {}
FESubplexOptimizeMethod::FESubplexOptimizeMethod(FEModel* fem) : NLOptOptimizeMethod(fem, NLOPT_LN_SBPLX) {}

#else
NLOptOptimizeMethod::NLOptOptimizeMethod(FEModel* fem, int algo) : FEOptimizeMethod(fem) {}
bool NLOptOptimizeMethod::Solve(FEOptimizeData* pOpt, vector<double>& amin, vector<double>& ymin, double* minObj) { return false; }
double NLOptOptimizeMethod::EvalObjective(unsigned int n, const double* x, double* grad) { return 0.0; }
void NLOptOptimizeMethod::EvalApproxGradient(unsigned int n, const double* x, double* grad) {}

FEMMAOptimizeMethod::FEMMAOptimizeMethod(FEModel* fem) : NLOptOptimizeMethod(fem, -1) {}
FENelderMeadOptMethod::FENelderMeadOptMethod(FEModel* fem) : NLOptOptimizeMethod(fem, -1) {}
FESubplexOptimizeMethod::FESubplexOptimizeMethod(FEModel* fem) : NLOptOptimizeMethod(fem, -1) {}

#endif
