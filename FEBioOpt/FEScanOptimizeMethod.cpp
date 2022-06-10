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
#include "FEScanOptimizeMethod.h"
#include "FEOptimizeData.h"
#include "FECore/log.h"

BEGIN_FECORE_CLASS(FEScanOptimizeMethod, FEOptimizeMethod)

END_FECORE_CLASS();

FEScanOptimizeMethod::FEScanOptimizeMethod(FEModel* fem) : FEOptimizeMethod(fem)
{

}

bool FEScanOptimizeMethod::Solve(FEOptimizeData* pOpt, vector<double>& amin, vector<double>& ymin, double* minObj)
{
	if (pOpt == 0) return false;
	FEOptimizeData& opt = *pOpt;
	FEObjectiveFunction& obj = opt.GetObjective();

	// set the intial values for the variables
	int ma = opt.InputParameters();
	vector<double> a(ma);
	for (int i=0; i<ma; ++i)
	{
		FEInputParameter* var = opt.GetInputParameter(i);
		a[i] = var->MinValue();
	}

	// loop until done
	vector<double> y(ma, 0.0);
	bool bdone = false;
	double fmin = 0.0;
	do
	{
		// solve the problem with the new input parameters
		if (opt.FESolve(a) == false) return false;

		// calculate objective function
		double fobj = obj.Evaluate(y);

		// update minimum
		if ((fmin == 0.0) || (fobj < fmin))
		{
			fmin = fobj;
			amin = a;
			ymin = y;
		}

		// update indices
		for (int i=0; i<ma; ++i)
		{
			FEInputParameter& vi = *opt.GetInputParameter(i);
			a[i] += vi.ScaleFactor();
			if (a[i] <= vi.MaxValue()) break;
			else if (i<ma-1) a[i] = vi.MinValue();
			else { bdone = true; }
		}
	}
	while (!bdone);

	// store the optimum data
	if (minObj) *minObj = fmin;

	return true;
}
