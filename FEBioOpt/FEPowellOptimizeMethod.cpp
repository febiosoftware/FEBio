#include "stdafx.h"
#include "FEPowellOptimizeMethod.h"
#include "FEOptimizeData.h"
#include <FECore/FEAnalysis.h>
#include "FECore/log.h"
#include <FECore/tools.h>

FEPowellOptimizeMethod* FEPowellOptimizeMethod::m_pThis = 0;

//-----------------------------------------------------------------------------
// FEPowellOptimizeMethod
//-----------------------------------------------------------------------------

bool FEPowellOptimizeMethod::Solve(FEOptimizeData *pOpt, vector<double>& amin, vector<double>& ymin, double* minObj)
{
	m_pOpt = pOpt;
	FEOptimizeData& opt = *pOpt;

	int nvar = opt.InputParameters();

	// set the initial guess
	vector<double> p(nvar);
	for (int i=0; i<nvar; ++i) p[i] = opt.GetInputParameter(i)->GetValue();

	// set the initial search directions
	vector<double> xi(nvar*nvar);
	for (int i=0; i<nvar; ++i)
	{
		for (int j=0; j<nvar; ++j) xi[i*nvar + j] = 0;
		xi[i*nvar + i] = 0.05*(1.0 + p[i]);
	}

	// don't forget to set this
	m_pThis = this;

	// call the powell routine
	int niter = 0;
	double fret = 0;
	powell(&p[0], &xi[0], nvar, 0.001, &niter, &fret, objfun);

	// store optimal values
	amin = p;
	if (minObj) *minObj = fret;

	return true;
}

//------------------------------------------------------------------

double FEPowellOptimizeMethod::ObjFun(double *p)
{
	// get the optimization data
	FEOptimizeData& opt = *m_pOpt;

	// set the input parameters
	int nvar = opt.InputParameters();
	vector<double> a(nvar);
	for (int i=0; i<nvar; ++i) a[i] = p[i];

	// solve the FE problem with the new parameters
	if (opt.FESolve(a) == false)
	{
		felog.printf("\n\n\nAAAAAAAAARRRRRRRRRGGGGGGGGHHHHHHHHHHH !!!!!!!!!!!!!\n\n\n\n");
		return 0;
	}
	else
	{
		// evaluate objective function
		FEObjectiveFunction& obj = opt.GetObjective();
		double fobj = obj.Evaluate();
		return fobj;
	}
}
