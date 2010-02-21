// LoadCurve.cpp: implementation of the LoadCurve class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "LoadCurve.h"

//-----------------------------------------------------------------------------
// FUNCTION : LoadCurve::Create
// Creates the time and data value arrays
//
void FELoadCurve::Create(int n)
{
	m_lp.resize(n);
}

//-----------------------------------------------------------------------------
// FUNCTION : LoadCurve::SetPoint
// Sets the time and data value of point i
// This function assumes that the load curve data has already been created
//
void FELoadCurve::SetPoint(int i, double time, double val)
{
	m_lp[i].time  = time;
	m_lp[i].value = val;
}

//-----------------------------------------------------------------------------
//! This function adds a datapoint to the loadcurve. The datapoint is inserted
//! at the appropriate place by examining the time parameter.

void FELoadCurve::Add(double time, double value)
{
	// find the place to insert the data point
	int n = 0;
	int nsize = Points();
	while ((n<nsize) && (m_lp[n].time < time)) ++n;

	// create a new data point
	LOADPOINT p = {time, value};

	// insert loadpoint
	m_lp.insert(p, n);
}

//-----------------------------------------------------------------------------
// FUNCTION : LoadCurve::Value
// Returns the load curve's value at time t.
// When the time value is outside the time range, the return value
// is that of the closest data value.
//
// TODO: maybe I should extrapolate the out-of-domain return values,
// in stead of clamping them. I think that is what NIKE does. Or even
// better let the user determine the out-of-range behaviour. Options could
// be zero, clamp to range, linear extrapolation, ...
//

double FELoadCurve::Value(double t)
{
	int n = Points();

	LOADPOINT* lp = m_lp;

	double value = 0;

	if (t < lp[0].time) value = lp[0].value;
	else if (t >= lp[n-1].time) value = lp[n-1].value;
	else
	{
		// find interval where t lies in
		// on break t lies in interval [l-1,l[
		int l=0;
		while (t >= lp[l].time) ++l;

		// calculate return value based on interpolation function
		switch (m_fnc)
		{
		case STEP:
			value = lp[l].value;
			break;
		case LINEAR:
			value = lp[l-1].value + (lp[l].value - lp[l-1].value)*(t - lp[l-1].time)/(lp[l].time - lp[l-1].time);
			break;
		}
	}

	return value;
}

//-----------------------------------------------------------------------------
// FUNCTION : LoadCurve::FindPoint(double t)
// This function finds the index of the first load point 
// for which the time is greater than t.
// It returns -1 if t is larger than the last time value
//

int FELoadCurve::FindPoint(double t)
{
	for (int i=0; i<Points(); ++i) if (m_lp[i].time > t) return i;
	return -1;
}

//-----------------------------------------------------------------------------

bool FELoadCurve::HasPoint(double t)
{
	const double tmax = m_lp[Points()-1].time;
	const double eps = 1e-7 * tmax;

	for (int i=0; i<Points(); ++i) if (fabs(m_lp[i].time - t) < eps) return true;

	return false;
}

//-----------------------------------------------------------------------------

void FELoadCurve::Serialize(Archive &ar)
{
	if (ar.IsSaving())
	{
		int n = Points();
		ar << n;
		for (int j=0; j<n; ++j)
		{
			LOADPOINT& p = LoadPoint(j);
			ar << p.time << p.value;
		}
	}
	else
	{
		int j, n;
		ar >> n;
		Create(n);
		for (j=0; j<n; ++j)
		{
			LOADPOINT& p = LoadPoint(j);
			ar >> p.time >> p.value;
		}
	}
}
