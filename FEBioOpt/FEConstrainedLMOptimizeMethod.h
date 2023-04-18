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



#pragma once
#include "FEOptimizeMethod.h"
#include <FECore/matrix.h>

//----------------------------------------------------------------------------
//! Optimization method using contrained Levenberg-Marquardt method
#ifdef HAVE_LEVMAR
class FEConstrainedLMOptimizeMethod : public FEOptimizeMethod
{
public:
	FEConstrainedLMOptimizeMethod(FEModel* fem);
	bool Solve(FEOptimizeData* pOpt, vector<double>& amin, vector<double>& ymin, double* minObj) override;

	FEOptimizeData* GetOptimizeData() { return m_pOpt; }

protected:
	FEOptimizeData* m_pOpt;

	void ObjFun(double* p, double* hx, int m, int n);

	static void objfun(double* p, double* hx, int m, int n, void* adata) 
	{ 
		FEConstrainedLMOptimizeMethod* clm = (FEConstrainedLMOptimizeMethod*)adata;
		return clm->ObjFun(p, hx, m , n);
	}

public:
	double	m_tau;		// scale factor for mu
	double	m_objtol;	// objective tolerance
	double	m_fdiff;	// forward difference step size
	int		m_nmax;		// maximum number of iterations
    int     m_loglevel; // log file output level
	bool	m_scaleParams;	// scale parameters flag

public:
	vector<double>	m_yopt;	// optimal y-values

	DECLARE_FECORE_CLASS();
};
#endif
