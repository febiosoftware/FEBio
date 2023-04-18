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
#include "FEParameterSweep.h"
#include <XML/XMLReader.h>
#include <FECore/FEModel.h>
#include <FECore/FEAnalysis.h>
#include <FECore/log.h>

FESweepParam::FESweepParam()
{
	m_min = m_max = m_step = 0.0;
	m_pd = nullptr;
}

FESweepParam::FESweepParam(const FESweepParam& p)
{
	m_min = p.m_min;
	m_max = p.m_max;
	m_step = p.m_step;
	m_paramName = p.m_paramName;
	m_pd = p.m_pd;
}

void FESweepParam::operator = (const FESweepParam& p)
{
	m_min = p.m_min;
	m_max = p.m_max;
	m_step = p.m_step;
	m_paramName = p.m_paramName;
	m_pd = p.m_pd;
}

void FESweepParam::SetValue(double v)
{
	*m_pd = v;
}

FEParameterSweep::FEParameterSweep(FEModel* fem) : FECoreTask(fem)
{
	m_niter = 0;
}

//! initialization
bool FEParameterSweep::Init(const char* szfile)
{
	// read the control file
	if (Input(szfile) == false) return false;

	// initialize the model
	if (GetFEModel()->Init() == false) return false;

	// check the parameters
	if (InitParams() == false) return false;

	// don't plot anything
	GetFEModel()->GetCurrentStep()->SetPlotHint(FE_PLOT_APPEND);
	GetFEModel()->GetCurrentStep()->SetPlotLevel(FE_PLOT_FINAL);

	return true;
}

bool FEParameterSweep::InitParams()
{
	// Make sure we have something to do
	if (m_params.empty()) return false;

	// check the parameters
	FEModel& fem = *GetFEModel();
	for (size_t i = 0; i<m_params.size(); ++i)
	{
		FESweepParam& p = m_params[i];

		// find the variable
		string name = p.m_paramName;
		FEParamValue val = fem.GetParameterValue(ParamString(name.c_str()));

		// see if we found the parameter
		if (val.isValid() == false)
		{
			feLogError("Cannot find parameter %s", name.c_str());
			return false;
		}

		// see if it's the correct type
		if (val.type() != FE_PARAM_DOUBLE)
		{
			feLogError("Invalid parameter type for parameter %s", name.c_str());
			return false;
		}

		// make sure we have a valid data pointer
		double* pd = (double*)val.data_ptr();
		if (pd == 0)
		{
			feLogError("Invalid data pointer for parameter %s", name.c_str());
			return false;
		}

		// store the pointer to the parameter
		p.m_pd = pd;
	}

	return true;
}

bool FEParameterSweep::Input(const char* szfile)
{
	// open the xml file
	XMLReader xml;
	if (xml.Open(szfile) == false) return false;

	// find the root tag
	XMLTag tag;
	if (xml.FindTag("febio_sweep", tag) == false) return false;

	++tag;
	do
	{
		if (tag == "param")
		{
			FESweepParam p;

			// read the parameter name
			const char* szname = tag.AttributeValue("name");
			p.m_paramName = szname;

			// read the values
			double d[3];
			int n = tag.value(d, 3);
			if (n != 3) throw XMLReader::InvalidValue(tag);

			// some error checking
			p.m_min = d[0];
			p.m_max = d[1];
			p.m_step = d[2];
			if (p.m_max < p.m_min) throw XMLReader::InvalidValue(tag);
			if (p.m_step <= 0.0) throw XMLReader::InvalidValue(tag);

			// looks good, so throw it on the pile
			m_params.push_back(p);
		}
		else throw XMLReader::InvalidTag(tag);
		++tag;
	} while (!tag.isend());

	// all done
	xml.Close();

	return true;
}

//! Run the optimization module
bool FEParameterSweep::Run()
{
	size_t ma = m_params.size();
	vector<double> a(ma);
	for (size_t i = 0; i<ma; ++i)
	{
		FESweepParam& pi = m_params[i];
		a[i] = pi.m_min;
	}

	// run the parameter sweep
	bool bdone = false;
	do
	{
		// solve the problem with the new input parameters
		if (FESolve(a) == false) return false;

		// update indices
		for (size_t i = 0; i<ma; ++i)
		{
			FESweepParam& pi = m_params[i];
			a[i] += pi.m_step;
			if (a[i] <= pi.m_max) break;
			else if (i<ma - 1) a[i] = pi.m_min;
			else { bdone = true; }
		}
	}
	while (!bdone);

	return true;
}

bool FEParameterSweep::FESolve(const vector<double>& a)
{
	++m_niter;
	feLog("\n----- Iteration: %d -----\n", m_niter);

	// set the input parameters
	size_t nvar = m_params.size();
	assert(nvar == a.size());
	for (int i = 0; i<nvar; ++i)
	{
		FESweepParam& var = m_params[i];
		var.SetValue(a[i]);

		string name = var.m_paramName;
		feLog("%-15s = %lg\n", name.c_str(), a[i]);
	}

	// reset the FEM data
	FEModel& fem = *GetFEModel();
	fem.BlockLog();

	// reset model
	fem.Reset();

	// solve the FE problem
	bool bret = fem.Solve();

	fem.UnBlockLog();

	return bret;
}
