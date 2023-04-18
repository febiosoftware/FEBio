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
#include "FEPowellOptimizeMethod.h"
#include "FEOptimizeData.h"
#include <FECore/FEAnalysis.h>
#include <FECore/log.h>
#include <FECore/tools.h>

FEPowellOptimizeMethod* FEPowellOptimizeMethod::m_pThis = 0;

FEPowellOptimizeMethod::FEPowellOptimizeMethod(FEModel* fem) : FEOptimizeMethod(fem)
{

}

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
	FEModel* fem = opt.GetFEModel();

	// set the input parameters
	int nvar = opt.InputParameters();
	vector<double> a(nvar);
	for (int i=0; i<nvar; ++i) a[i] = p[i];

	// solve the FE problem with the new parameters
	if (opt.FESolve(a) == false)
	{
		feLogEx(fem, "\n\n\nAAAAAAAAARRRRRRRRRGGGGGGGGHHHHHHHHHHH !!!!!!!!!!!!!\n\n\n\n");
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
