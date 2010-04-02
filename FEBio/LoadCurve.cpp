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


inline double lerp(double t, double t0, double f0, double t1, double f1)
{
	return f0 + (f1 - f0)*(t - t0)/(t1 - t0);
}

inline double qerp(double t, double t0, double f0, double t1, double f1, double t2, double f2)
{
	double q0 = ((t2 - t )*(t1 - t ))/((t2 - t0)*(t1 - t0));
	double q1 = ((t2 - t )*(t  - t0))/((t2 - t1)*(t1 - t0));
	double q2 = ((t  - t1)*(t  - t0))/((t2 - t1)*(t2 - t0));

	return f0*q0 + f1*q1 + f2*q2;
}

double FELoadCurve::Value(double time)
{
	LOADPOINT* lp = m_lp;

	int nsize = Points();
	if (nsize == 0) return 0;
	if (nsize == 1) return lp[0].value;

	int N = nsize - 1;

	if (time == lp[0].time) return lp[0].value;
	if (time == lp[N].time) return lp[N].value;

	if (time < lp[0].time) return ExtendValue(time);
	if (time > lp[N].time) return ExtendValue(time);


	if (m_fnc == LINEAR)
	{
		int n = 0;
		while (lp[n].time <= time) ++n;
	
		double t0 = lp[n-1].time;
		double t1 = lp[n  ].time;

		double f0 = lp[n-1].value;
		double f1 = lp[n  ].value;

		return lerp(time, t0, f0, t1, f1);
	}
	else if (m_fnc == STEP)
	{
		int n=0;
		while (lp[n].time < time) ++n;

		return lp[n].value;
	}
	else if (m_fnc == SMOOTH)
	{
		if (nsize == 2)
		{
			double t0 = lp[0].time;
			double t1 = lp[1].time;

			double f0 = lp[0].value;
			double f1 = lp[1].value;

			return lerp(time, t0, f0, t1, f1);
		}
		else if (nsize == 3)
		{
			double t0 = lp[0].time;
			double t1 = lp[1].time;
			double t2 = lp[2].time;

			double f0 = lp[0].value;
			double f1 = lp[1].value;
			double f2 = lp[2].value;

			return qerp(time, t0, f0, t1, f1, t2, f2);
		}
		else
		{
			int n = 0;
			while (lp[n].time <= time) ++n;

			if (n == 1)
			{
				double t0 = lp[0].time;
				double t1 = lp[1].time;
				double t2 = lp[2].time;

				double f0 = lp[0].value;
				double f1 = lp[1].value;
				double f2 = lp[2].value;

				return qerp(time, t0, f0, t1, f1, t2, f2);
			}
			else if (n == nsize-1)
			{
				double t0 = lp[n-2].time;
				double t1 = lp[n-1].time;
				double t2 = lp[n  ].time;

				double f0 = lp[n-2].value;
				double f1 = lp[n-1].value;
				double f2 = lp[n  ].value;

				return qerp(time, t0, f0, t1, f1, t2, f2);
			}
			else
			{
				double t0 = lp[n-2].time;
				double t1 = lp[n-1].time;
				double t2 = lp[n  ].time;
				double t3 = lp[n+1].time;

				double f0 = lp[n-2].value;
				double f1 = lp[n-1].value;
				double f2 = lp[n  ].value;
				double f3 = lp[n+1].value;

				double q1 = qerp(time, t0, f0, t1, f1, t2, f2);
				double q2 = qerp(time, t1, f1, t2, f2, t3, f3);

				return lerp(time, t1, q1, t2, q2);
			}
		}
	}

	return 0;
}

//-----------------------------------------------------------------------------
//! This function determines the value of the load curve outside of its domain
//!
double FELoadCurve::ExtendValue(double t)
{
	LOADPOINT* lp = m_lp;

	int nsize =Points();
	int N = nsize - 1;

	if (nsize == 0) return 0;
	if (nsize == 1) return lp[0].value;

	double Dt = (lp[N].time - lp[0].time);
	double dt = 0.001*Dt;
	if (dt == 0) return lp[0].value;

	switch (m_ext)
	{
	case CONSTANT:
		if (t < lp[0].time) return lp[0].value;
		if (t > lp[N].time) return lp[N].value;
		break;
	case EXTRAPOLATE:
		switch (m_fnc)
		{
		case STEP:
			{
				if (t < lp[0].time) return lp[0].value;
				if (t > lp[N].time) return lp[N].value;
			}
			break;
		case LINEAR:
			{
				if (t < lp[0].time) return lerp(t, lp[0].time, lp[0].value, lp[1].time, lp[1].value);
				else return lerp(t, lp[N-1].time, lp[N-1].value, lp[N].time, lp[N].value);
			}
			break;
		case SMOOTH:
			{
				if (t < lp[0].time) return lerp(t, lp[0].time, lp[0].value, lp[0].time + dt, Value(lp[0].time+dt));
				else return lerp(t, lp[N].time - dt, Value(lp[N].time - dt), lp[N].time, lp[N].value);
			}
			return 0;
		}
		break;
	case REPEAT:
		{
			if (t < lp[0].time) while (t < lp[0].time) t += Dt;
			else while (t > lp[N].time) t -= Dt;
			return Value(t);
		}
		break;
	case REPEAT_OFFSET:
		{
			int n = 0;
			if (t < lp[0].time) while (t < lp[0].time) { t += Dt; --n; }
			else while (t > lp[N].time) { t -= Dt; ++n; }
			double off = n*(lp[N].value - lp[0].value);

			return Value(t)+off;
		}
		break;
	}

	return 0;
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
