#include "stdafx.h"
#include "FEObjectiveFunction.h"
#include <FECore/FEModel.h>
#include <FECore/log.h>

//=============================================================================

FEObjectiveFunction::FEObjectiveFunction(FEModel* fem) : m_fem(fem)
{
	m_verbose = true;
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

	double chisq = 0.0;
	if (m_verbose) felog.printf("               CURRENT        REQUIRED      DIFFERENCE\n");
	for (int i = 0; i<ndata; ++i)
	{
		double dy = (y[i] - y0[i]);
		chisq += dy*dy;
		if (m_verbose) felog.printf("%5d: %15.10lg %15.10lg %15lg\n", i + 1, y[i], y0[i], fabs(y[i] - y0[i]));
	}
	felog.printf("objective value: %lg\n", chisq);

	return chisq;
}

//=============================================================================

//----------------------------------------------------------------------------
FEDataFitObjective::FEDataFitObjective(FEModel* fem) : FEObjectiveFunction(fem), m_lc(fem)
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
	FEModel& fem = *GetFEM();

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
	for (int i = 0; i<ndata; ++i) y0[i] = m_lc.LoadPoint(i).value;
}

//----------------------------------------------------------------------------
void FEDataFitObjective::EvaluateFunctions(vector<double>& f)
{
	int ndata = m_lc.Points();
	for (int i = 0; i<ndata; ++i)
	{
		double ti = m_lc.LoadPoint(i).time;
		f[i] = m_src->Evaluate(ti);
	}
}

//=============================================================================
FEMinimizeObjective::FEMinimizeObjective(FEModel* fem) : FEObjectiveFunction(fem)
{
}

bool FEMinimizeObjective::AddFunction(const char* szname, double targetValue)
{
	// make sure we have a model
	FEModel* fem = GetFEM();
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

	FEModel* fem = GetFEM();
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
