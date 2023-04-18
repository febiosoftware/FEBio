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
#include "FEObjectiveFunction.h"
#include <FECore/FEModel.h>
#include <FECore/log.h>

//=============================================================================

FEObjectiveFunction::FEObjectiveFunction(FEModel* fem) : m_fem(fem)
{
	m_verbose = false;
}

FEObjectiveFunction::~FEObjectiveFunction()
{
	
}

bool FEObjectiveFunction::Init()
{
	// make sure we have a model
	if (m_fem == 0) return false;

	return true;
}

void FEObjectiveFunction::Reset()
{
}

double FEObjectiveFunction::Evaluate()
{
	vector<double> dummy(Measurements());
	return Evaluate(dummy);
}

double FEObjectiveFunction::Evaluate(vector<double>& y)
{
	// get the number of measurements
	int ndata = Measurements();
	y.resize(ndata);

	// evaluate the functions
	EvaluateFunctions(y);

	// get the measurement vector
	vector<double> y0(ndata);
	GetMeasurements(y0);
    
    // evaluate regression coefficient R^2
	double rsq = RegressionCoefficient(y0, y);

	double chisq = 0.0;
	if (m_verbose) feLog("               CURRENT        REQUIRED      DIFFERENCE\n");
	for (int i = 0; i<ndata; ++i)
	{
		double dy = (y[i] - y0[i]);
		chisq += dy*dy;
		if (m_verbose) feLog("%5d: %15.10lg %15.10lg %15lg\n", i + 1, y[i], y0[i], fabs(y[i] - y0[i]));
	}
	feLog("objective value: %lg\n", chisq);
    feLog("regression coef: %lg\n", rsq);

	return chisq;
}

double FEObjectiveFunction::RegressionCoefficient(const std::vector<double>& y0, const std::vector<double>& y)
{
	int ndata = (int)y0.size();
	double xb = 0, yb = 0, xyb = 0, x2b = 0, y2b = 0;
	for (int i = 0; i < ndata; ++i)
	{
		xb += y0[i]; yb += y[i];
		xyb += y0[i] * y[i];
		x2b += pow(y0[i], 2); y2b += pow(y[i], 2);
	}
	xb /= ndata; yb /= ndata;
	xyb /= ndata;
	x2b /= ndata; y2b /= ndata;

	double rsq = pow(xyb - xb * yb, 2) / (x2b - xb * xb) / (y2b - yb * yb);
	return rsq;
}

//=============================================================================

//----------------------------------------------------------------------------
FEDataFitObjective::FEDataFitObjective(FEModel* fem) : FEObjectiveFunction(fem)
{
	m_src = 0;
}

FEDataFitObjective::~FEDataFitObjective()
{
	if (m_src) delete m_src;
	m_src = 0;
}

//----------------------------------------------------------------------------
bool FEDataFitObjective::Init()
{
	if (FEObjectiveFunction::Init() == false) return false;

	// get the FE model
	FEModel& fem = *GetFEModel();

	// initialize data source
	if (m_src == 0) return false;
	if (m_src->Init() == false) return false;

	return true;
}

//----------------------------------------------------------------------------
// set the data source
void FEDataFitObjective::SetDataSource(FEDataSource* src)
{
	if (m_src) delete m_src;
	m_src = src;
}

//----------------------------------------------------------------------------
// set the data measurements
void FEDataFitObjective::SetMeasurements(const vector<pair<double, double> >& data)
{
	m_lc.Clear();
	int n = (int)data.size();
	for (int i=0; i<n; ++i)
	{
		const pair<double,double>& pt = data[i];
		m_lc.Add(pt.first, pt.second);
	}
}

//----------------------------------------------------------------------------
void FEDataFitObjective::Reset()
{
	// call base class first
	FEObjectiveFunction::Reset();

	m_src->Reset();
}

//----------------------------------------------------------------------------
// return the number of measurements. I.e. the size of the measurement vector
int FEDataFitObjective::Measurements()
{
	return m_lc.Points();
}

//----------------------------------------------------------------------------
// Evaluate the measurement vector and return in y0
void FEDataFitObjective::GetMeasurements(vector<double>& y0)
{
	int ndata = m_lc.Points();
	y0.resize(ndata);
	for (int i = 0; i<ndata; ++i) y0[i] = m_lc.Point(i).y();
}

//----------------------------------------------------------------------------
void FEDataFitObjective::EvaluateFunctions(vector<double>& f)
{
	int ndata = m_lc.Points();
	for (int i = 0; i<ndata; ++i)
	{
		double xi = m_lc.Point(i).x();
		f[i] = m_src->Evaluate(xi);
	}
}

//=============================================================================
FEMinimizeObjective::FEMinimizeObjective(FEModel* fem) : FEObjectiveFunction(fem)
{
}

bool FEMinimizeObjective::AddFunction(const char* szname, double targetValue)
{
	// make sure we have a model
	FEModel* fem = GetFEModel();
	if (fem == 0) return false;

	// create a new function
	Function fnc;
	fnc.name = szname;
	fnc.y0 = targetValue;

	m_Func.push_back(fnc);

	return true;
}

bool FEMinimizeObjective::Init()
{
	if (FEObjectiveFunction::Init() == false) return false;

	FEModel* fem = GetFEModel();
	if (fem == 0) return false;

	int N = (int) m_Func.size();
	for (int i=0; i<N; ++i)
	{
		Function& Fi = m_Func[i];
		FEParamValue val = fem->GetParameterValue(ParamString(Fi.name.c_str()));
		if (val.isValid() == false) return false;
		if (val.type() != FE_PARAM_DOUBLE) return false;
		Fi.var = (double*) val.data_ptr();
		if (Fi.var == 0) return false;
	}

	return true;
}

int FEMinimizeObjective::Measurements()
{
	return (int) m_Func.size();
}


void FEMinimizeObjective::EvaluateFunctions(vector<double>& f)
{
	int N = (int)m_Func.size();
	f.resize(N);
	for (int i=0; i<N; ++i)
	{
		Function& Fi = m_Func[i];
		if (Fi.var) f[i] = *Fi.var;
		else f[i] = 0;
	}
}

void FEMinimizeObjective::GetMeasurements(vector<double>& y)
{
	int N = (int) m_Func.size();
	y.resize(N);
	for (int i=0; i<N; ++i)
	{
		y[i] = m_Func[i].y0;
	}
}

//=============================================================================
FEElementDataTable::FEElementDataTable(FEModel* fem) : FEObjectiveFunction(fem)
{
	m_var = nullptr;
}

void FEElementDataTable::AddValue(int elemID, double v)
{
	Entry d;
	d.elemId = elemID;
	d.target = v;
	d.pe = nullptr;
	m_Data.push_back(d);
}

void FEElementDataTable::SetVariable(FELogElemData* var)
{
	m_var = var;
}

bool FEElementDataTable::Init()
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	int N = (int)m_Data.size();
	if (N == 0) return false;

	for (int i = 0; i < N; ++i)
	{
		Entry& di = m_Data[i];
		FEElement* el = mesh.FindElementFromID(di.elemId);
		if (el == nullptr) return false;
		di.pe = el;
	}

	return true;
}

// return number of measurements (i.e. nr of terms in objective function)
int FEElementDataTable::Measurements()
{
	return (int)m_Data.size();
}

// evaluate the function values (i.e. the f_i above)
void FEElementDataTable::EvaluateFunctions(vector<double>& f)
{
	assert(m_var);

	int N = (int)m_Data.size();
	f.resize(N);
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	for (int i = 0; i < N; ++i)
	{
		FEElement* pe = m_Data[i].pe;

		// calculate element average measure
		double val = m_var->value(*pe);

		// store result
		f[i] = val;
	}
}

// get the measurement vector (i.e. the y_i above)
void FEElementDataTable::GetMeasurements(vector<double>& y)
{
	int N = (int)m_Data.size();
	y.resize(N);
	for (int i = 0; i < N; ++i)
	{
		y[i] = m_Data[i].target;
	}
}


//=============================================================================
FENodeDataTable::FENodeDataTable(FEModel* fem) : FEObjectiveFunction(fem)
{
	
}

bool FENodeDataTable::AddValue(int nodeID, vector<double>& v)
{
	if (v.size() != m_var.size()) return false;

	for (int i = 0; i < v.size(); ++i)
	{
		Entry d;
		d.nodeId = nodeID;
		d.target = v[i];
		d.ivar = i;
		d.index = -1;
		m_Data.push_back(d);
	}

	return true;
}

void FENodeDataTable::AddVariable(FELogNodeData* var)
{
	m_var.push_back(var);
}

bool FENodeDataTable::Init()
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	int N = (int)m_Data.size();
	if (N == 0) return false;

	for (int i = 0; i < N; ++i)
	{
		Entry& di = m_Data[i];
		FENode* node = mesh.FindNodeFromID(di.nodeId);

		if (node == nullptr) return false;
		di.index = di.nodeId - 1;	// NOTE: This assumes one-based indexing of node IDs. 
	}

	return true;
}

// return number of measurements (i.e. nr of terms in objective function)
int FENodeDataTable::Measurements()
{
	return (int)m_Data.size();
}

// evaluate the function values (i.e. the f_i above)
void FENodeDataTable::EvaluateFunctions(vector<double>& f)
{
	assert(m_var.size() > 0);

	int N = (int)m_Data.size();
	f.resize(N);
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	for (int i = 0; i < N; ++i)
	{
		int n = m_Data[i].index;
		int v = m_Data[i].ivar;

		// calculate node value
		double val = m_var[v]->value(mesh.Node(n));

		// store result
		f[i] = val;
	}
}

// get the measurement vector (i.e. the y_i above)
void FENodeDataTable::GetMeasurements(vector<double>& y)
{
	int N = (int)m_Data.size();
	y.resize(N);
	for (int i = 0; i < N; ++i)
	{
		y[i] = m_Data[i].target;
	}
}
