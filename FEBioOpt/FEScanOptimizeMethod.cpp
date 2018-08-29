#include "stdafx.h"
#include "FEScanOptimizeMethod.h"
#include "FEOptimizeData.h"
#include "FECore/log.h"

BEGIN_PARAMETER_LIST(FEScanOptimizeMethod, FEOptimizeMethod)

END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
// FEScanOptimizeMethod
//-----------------------------------------------------------------------------

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
