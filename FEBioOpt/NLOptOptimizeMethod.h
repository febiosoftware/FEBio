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
#pragma once
#include "FEOptimizeMethod.h"
#include "febioopt_api.h"

// Base class for optimization methods that use the NLOpt library.
// This class provides a common implementation of the Solve function that sets up the optimization problem and calls the NLOpt solver. Derived classes can specify the particular NLOpt algorithm to use and can implement the EvalObjective function to evaluate the objective function and its gradient.
class FEBIOOPT_API NLOptOptimizeMethod : public FEOptimizeMethod
{
public:
	NLOptOptimizeMethod(FEModel* fem, int nlopt_algorithm);
	bool Solve(FEOptimizeData* pOpt, vector<double>& amin, vector<double>& ymin, double* minObj) override;

	double EvalObjective(unsigned int n, const double* x, double* grad);

private:
	void EvalApproxGradient(unsigned int n, const double* x, double* grad);

private:
	int m_nlopt_algorithm = -1;
	FEOptimizeData* m_optData = nullptr;
	int m_max_iter = 1000;
	double m_ftol_abs = 0; // use nlopt's default ftol_abs
	double m_ftol_rel = 0; // use nlopt's default ftol_rel

	DECLARE_FECORE_CLASS();
};

class FEBIOOPT_API FEMMAOptimizeMethod : public NLOptOptimizeMethod
{
public:
	FEMMAOptimizeMethod(FEModel* fem);
};

class FEBIOOPT_API FENelderMeadOptMethod : public NLOptOptimizeMethod
{
public:
	FENelderMeadOptMethod(FEModel* fem);
};

class FEBIOOPT_API FESubplexOptimizeMethod : public NLOptOptimizeMethod
{
public:
	FESubplexOptimizeMethod(FEModel* fem);
};
