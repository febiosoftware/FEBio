#include "stdafx.h"
#include "FEObjectiveFunction.h"
#include <FECore/FEModel.h>
#include <FECore/log.h>
#include <FEBioMech/FEElasticMaterial.h>

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
		double xi = m_lc.LoadPoint(i).time;
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

//=============================================================================
FEElementDataTable::FEElementDataTable(FEModel* fem) : FEObjectiveFunction(fem)
{
	m_var = 0;
}

void FEElementDataTable::AddValue(int elemID, double v)
{
	Entry d;
	d.elemId = elemID;
	d.target = v;
	d.pe = nullptr;
	m_Data.push_back(d);
}

void FEElementDataTable::SetVariable(int n)
{
	m_var = n;
}

bool FEElementDataTable::Init()
{
	FEModel& fem = *GetFEM();
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
	int N = (int)m_Data.size();
	f.resize(N);
	FEModel& fem = *GetFEM();
	FEMesh& mesh = fem.GetMesh();
	for (int i = 0; i < N; ++i)
	{
		FEElement* pe = m_Data[i].pe;

		double val = 0.0;
		
		// calculate element average measure
		if (m_var == 0)
		{
			int nint = pe->GaussPoints();
			mat3ds Eavg; Eavg.zero();
			for (int n = 0; n < nint; ++n)
			{
				FEMaterialPoint& mp = *pe->GetMaterialPoint(n);
				FEElasticMaterialPoint& ep = *mp.ExtractData<FEElasticMaterialPoint>();

				mat3ds C = ep.LeftCauchyGreen();
				mat3dd I(1.0);
				mat3ds E = (C - I)*0.5;

				Eavg += E;
			}
			Eavg /= (double)nint;
			val = Eavg.effective_norm();
		}
		else
		{
			int nint = pe->GaussPoints();
			mat3ds savg; savg.zero();
			for (int n = 0; n < nint; ++n)
			{
				FEMaterialPoint& mp = *pe->GetMaterialPoint(n);
				FEElasticMaterialPoint& ep = *mp.ExtractData<FEElasticMaterialPoint>();

				savg += ep.m_s;
			}
			savg /= (double)nint;
			val = savg.effective_norm();
		}

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
