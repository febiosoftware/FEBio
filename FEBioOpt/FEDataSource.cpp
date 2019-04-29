/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, The Trustees of Columbia University in 
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
#include "FEDataSource.h"
#include <FECore/FEModel.h>
#include <FECore/log.h>

//=================================================================================================
FEDataSource::FEDataSource(FEModel* fem) : m_fem(*fem)
{
	
}

FEDataSource::~FEDataSource()
{
	
}

bool FEDataSource::Init()
{
	return true;
}

void FEDataSource::Reset()
{
	
}

//=================================================================================================
bool FEDataParameter::update(FEModel* pmdl, unsigned int nwhen, void* pd)
{
	// get the optimizaton data
	FEDataParameter& src = *((FEDataParameter*)pd);
	src.update();

	return true;
}

void FEDataParameter::update()
{
	// get the current time value
	double time = m_fem.GetTime().currentTime;

	// evaluate the current reaction force value
	double x = *(m_px);
	double y = *(m_py);

	// add the data pair to the loadcurve
	m_rf.Add(x, y);
}

FEDataParameter::FEDataParameter(FEModel* fem) : FEDataSource(fem), m_rf(fem)
{
	m_ord = "fem.time";
}

void FEDataParameter::SetParameterName(const std::string& name)
{
	m_param = name;
}

// set the ordinate name
void FEDataParameter::SetOrdinateName(const std::string& name)
{
	m_ord = name;
}

bool FEDataParameter::Init()
{
	FEModel* fem = &m_fem;

	// find all the parameters
	FEParamValue val = m_fem.GetParameterValue(ParamString(m_param.c_str()));
	if (val.isValid() == false) {
		feLogErrorEx(fem, "Invalid parameter name %s", m_param.c_str());
		return false;
	}
	if (val.type() != FE_PARAM_DOUBLE) {
		feLogErrorEx(fem, "Invalid type for parameter %s", m_param.c_str());
		return false;
	}
	m_py = (double*)val.data_ptr();
	if (m_py == 0) {
		feLogErrorEx(fem, "Invalid data pointer for parameter %s", m_param.c_str());
		return false;
	}

	// find the ordinate
	val = m_fem.GetParameterValue(ParamString(m_ord.c_str()));
	if (val.isValid() == false) {
		feLogErrorEx(fem, "Invalid ordinate name %s", m_ord.c_str());
		return false;
	}
	if (val.type() != FE_PARAM_DOUBLE) {
		feLogErrorEx(fem, "Invalid type for ordinate %s", m_ord.c_str());
		return false;
	}
	m_px = (double*)val.data_ptr();
	if (m_px == 0) {
		feLogErrorEx(fem, "Invalid data pointer for ordinate %s", m_ord.c_str());
		return false;
	}

	// register callback
	m_fem.AddCallback(update, CB_INIT | CB_MAJOR_ITERS, (void*) this);

	return FEDataSource::Init();
}

void FEDataParameter::Reset()
{
	// reset the reaction force load curve
	m_rf.Clear();
	*m_px = 0.0;
	FEDataSource::Reset();
}

double FEDataParameter::Evaluate(double x)
{
	return m_rf.value(x);
}

//=================================================================================================
FEDataFilterPositive::FEDataFilterPositive(FEModel* fem) : FEDataSource(fem)
{
	m_src = 0;
}

FEDataFilterPositive::~FEDataFilterPositive()
{
	if (m_src) delete m_src;
}

void FEDataFilterPositive::SetDataSource(FEDataSource* src)
{
	if (m_src) delete m_src;
	m_src = src;
}

bool FEDataFilterPositive::Init()
{
	if (m_src == 0) return false;

	return m_src->Init();
}

void FEDataFilterPositive::Reset()
{
	if (m_src) m_src->Reset();
}

double FEDataFilterPositive::Evaluate(double t)
{
	double v = m_src->Evaluate(t);
	return (v >= 0.0 ? v : -v);
}
