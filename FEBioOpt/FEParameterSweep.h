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
#include <FECore/FECoreTask.h>

// This class represents a parameter that will be swept
class FESweepParam
{
public:
	FESweepParam();
	FESweepParam(const FESweepParam& p);
	void operator = (const FESweepParam& p);

	void SetValue(double v);

public:
	string	m_paramName;
	double	m_min, m_max, m_step;
	double*	m_pd;
};

// This task implements a parameter sweep, where the same model is run similar times,
// each time with one or more parameters modified.
class FEParameterSweep : public FECoreTask
{
public:
	FEParameterSweep(FEModel* fem);

	//! initialization
	bool Init(const char* szfile) override;

	//! Run the optimization module
	bool Run() override;

private:
	bool Input(const char* szfile);
	bool InitParams();
	bool FESolve(const vector<double>& a);

private:
	vector<FESweepParam>	m_params;
	int						m_niter;
};
