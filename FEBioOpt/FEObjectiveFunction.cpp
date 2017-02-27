#include "stdafx.h"
#include "FEObjectiveFunction.h"
#include <FECore/FEModel.h>
#include <FECore/log.h>

//=============================================================================

FEObjectiveFunction::FEObjectiveFunction(FEModel* fem) : m_fem(fem)
{

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
	felog.printf("               CURRENT        REQUIRED      DIFFERENCE\n");
	for (int i = 0; i<ndata; ++i)
	{
		double dy = (y[i] - y0[i]);
		chisq += dy*dy;
		felog.printf("%5d: %15.10lg %15.10lg %15lg\n", i + 1, y[i], y0[i], fabs(y[i] - y0[i]));
	}
	felog.printf("objective value: %lg\n", chisq);

	return chisq;
}

//=============================================================================

//-------------------------------------------------------------------
bool fecb(FEModel* pmdl, unsigned int nwhen, void* pd)
{
	// get the optimizaton data
	FEDataFitObjective& obj = *((FEDataFitObjective*)pd);

	// get the FEM data
	FEModel& fem = *obj.GetFEM();

	// get the current time value
	double time = fem.m_ftime;

	// evaluate the current reaction force value
	double value = *(obj.m_pd);

	// add the data pair to the loadcurve
	FELoadCurve& lc = obj.ReactionLoad();
	lc.Add(time, value);

	return true;
}

//----------------------------------------------------------------------------
FEDataFitObjective::FEDataFitObjective(FEModel* fem) : FEObjectiveFunction(fem)
{

}

//----------------------------------------------------------------------------
bool FEDataFitObjective::Init()
{
	if (FEObjectiveFunction::Init() == false) return false;

	// get the FE model
	FEModel& fem = *GetFEM();

	// find all the parameters
	m_pd = fem.FindParameter(m_szname);
	if (m_pd == 0) return false;

	// register callback
	fem.AddCallback(fecb, CB_INIT | CB_MAJOR_ITERS, (void*) this);

	return true;
}

//----------------------------------------------------------------------------
void FEDataFitObjective::Reset()
{
	// call base class first
	FEObjectiveFunction::Reset();

	// reset the reaction force load curve
	FELoadCurve& lc = ReactionLoad();
	lc.Clear();
}

//----------------------------------------------------------------------------
// return the number of measurements. I.e. the size of the measurement vector
int FEDataFitObjective::Measurements()
{
	FELoadCurve& lc = GetLoadCurve(m_nlc);
	return lc.Points();
}

//----------------------------------------------------------------------------
// Evaluate the measurement vector and return in y0
void FEDataFitObjective::GetMeasurements(vector<double>& y0)
{
	FELoadCurve& lc = GetLoadCurve(m_nlc);
	int ndata = lc.Points();
	y0.resize(ndata);
	for (int i = 0; i<ndata; ++i) y0[i] = lc.LoadPoint(i).value;
}

void FEDataFitObjective::EvaluateFunctions(vector<double>& f)
{
	FELoadCurve& lc = GetLoadCurve(m_nlc);
	int ndata = lc.Points();
	FELoadCurve& rlc = ReactionLoad();
	for (int i = 0; i<ndata; ++i)
	{
		double xi = lc.LoadPoint(i).time;
		f[i] = rlc.Value(xi);
	}
}
