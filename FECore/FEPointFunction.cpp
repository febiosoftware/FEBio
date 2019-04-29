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
#include "FEPointFunction.h"
#include "DumpStream.h"

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEPointFunction, FEFunction1D)
	ADD_PARAMETER(m_points, "point");
	ADD_PARAMETER(m_fnc, "interpolate", FE_PARAM_ATTRIBUTE, "step\0linear\0smooth\0");
	ADD_PARAMETER(m_ext, "extend"     , FE_PARAM_ATTRIBUTE, "constant\0extrapolate\0repeat\0repeat offset\0");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! default constructor
FEPointFunction::FEPointFunction(FEModel* fem) : FEFunction1D(fem), m_fnc(LINEAR), m_ext(CONSTANT)
{
}

//-----------------------------------------------------------------------------
//! Clears the loadcurve data
void FEPointFunction::Clear()
{ 
	m_points.clear();
}

//-----------------------------------------------------------------------------
//! return nr of points
int FEPointFunction::Points() const
{ 
	return (int) m_points.size(); 
}

//-----------------------------------------------------------------------------
// Sets the time and data value of point i
// This function assumes that the load curve data has already been created
//
void FEPointFunction::SetPoint(int i, double x, double y)
{
	vec2d& pt = m_points[i];
	pt.x() = x;
	pt.y() = y;
}

//-----------------------------------------------------------------------------
//! returns point i
LOADPOINT FEPointFunction::LoadPoint(int i) const
{ 
	const vec2d& p = m_points[i];
	LOADPOINT lp;
	lp.time  = p.x();
	lp.value = p.y();
	return lp; 
}

//-----------------------------------------------------------------------------
//! This function adds a datapoint to the loadcurve. The datapoint is inserted
//! at the appropriate place by examining the time parameter.

void FEPointFunction::Add(double x, double y)
{
	// find the place to insert the data point
	int n = 0;
	int nsize = Points();
	while ((n<nsize) && (m_points[n].x() < x)) ++n;

	// insert loadpoint
	m_points.insert(m_points.begin() + n, vec2d(x, y));
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
	return f0 + (f1 - f0)*(t - t0) / (t1 - t0);
}

inline double qerp(double t, double t0, double f0, double t1, double f1, double t2, double f2)
{
	double q0 = ((t2 - t)*(t1 - t)) / ((t2 - t0)*(t1 - t0));
	double q1 = ((t2 - t)*(t - t0)) / ((t2 - t1)*(t1 - t0));
	double q2 = ((t - t1)*(t - t0)) / ((t2 - t1)*(t2 - t0));

	return f0*q0 + f1*q1 + f2*q2;
}

double FEPointFunction::value(double time) const
{
	int nsize = Points();
	if (nsize == 0) return 0;
	if (nsize == 1) return m_points[0].y();

	int N = nsize - 1;

	if (time == m_points[0].x()) return m_points[0].y();
	if (time == m_points[N].x()) return m_points[N].y();

	if (time < m_points[0].x()) return ExtendValue(time);
	if (time > m_points[N].x()) return ExtendValue(time);

	if (m_fnc == LINEAR)
	{
		int n = 0;
		while (m_points[n].x() <= time) ++n;

		double t0 = m_points[n - 1].x();
		double t1 = m_points[n    ].x();

		double f0 = m_points[n - 1].y();
		double f1 = m_points[n    ].y();

		return lerp(time, t0, f0, t1, f1);
	}
	else if (m_fnc == STEP)
	{
		int n = 0;
		while (m_points[n].x() <= time) ++n;

		return m_points[n].y();
	}
	else if (m_fnc == SMOOTH)
	{
		if (nsize == 2)
		{
			double t0 = m_points[0].x();
			double t1 = m_points[1].x();

			double f0 = m_points[0].y();
			double f1 = m_points[1].y();

			return lerp(time, t0, f0, t1, f1);
		}
		else if (nsize == 3)
		{
			double t0 = m_points[0].x();
			double t1 = m_points[1].x();
			double t2 = m_points[2].x();

			double f0 = m_points[0].y();
			double f1 = m_points[1].y();
			double f2 = m_points[2].y();

			return qerp(time, t0, f0, t1, f1, t2, f2);
		}
		else
		{
			int n = 0;
			while (m_points[n].x() <= time) ++n;

			if (n == 1)
			{
				double t0 = m_points[0].x();
				double t1 = m_points[1].x();
				double t2 = m_points[2].x();

				double f0 = m_points[0].y();
				double f1 = m_points[1].y();
				double f2 = m_points[2].y();

				return qerp(time, t0, f0, t1, f1, t2, f2);
			}
			else if (n == nsize - 1)
			{
				double t0 = m_points[n - 2].x();
				double t1 = m_points[n - 1].x();
				double t2 = m_points[n    ].x();

				double f0 = m_points[n - 2].y();
				double f1 = m_points[n - 1].y();
				double f2 = m_points[n    ].y();

				return qerp(time, t0, f0, t1, f1, t2, f2);
			}
			else
			{
				double t0 = m_points[n - 2].x();
				double t1 = m_points[n - 1].x();
				double t2 = m_points[n    ].x();
				double t3 = m_points[n + 1].x();

				double f0 = m_points[n - 2].y();
				double f1 = m_points[n - 1].y();
				double f2 = m_points[n    ].y();
				double f3 = m_points[n + 1].y();

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
double FEPointFunction::ExtendValue(double t) const
{
	int nsize = Points();
	int N = nsize - 1;

	if (nsize == 0) return 0;
	if (nsize == 1) return m_points[0].y();

	double Dt = (m_points[N].x() - m_points[0].x());
	double dt = 0.001*Dt;
	if (dt == 0) return m_points[0].y();

	switch (m_ext)
	{
	case CONSTANT:
		if (t < m_points[0].x()) return m_points[0].y();
		if (t > m_points[N].x()) return m_points[N].y();
		break;
	case EXTRAPOLATE:
		switch (m_fnc)
		{
		case STEP:
		{
			if (t < m_points[0].x()) return m_points[0].y();
			if (t > m_points[N].x()) return m_points[N].y();
		}
			break;
		case LINEAR:
		{
			if (t < m_points[0].x()) return lerp(t, m_points[0].x(), m_points[0].y(), m_points[1].x(), m_points[1].y());
			else return lerp(t, m_points[N - 1].x(), m_points[N - 1].y(), m_points[N].x(), m_points[N].y());
		}
			break;
		case SMOOTH:
		{
			if (t < m_points[0].x()) return lerp(t, m_points[0].x(), m_points[0].y(), m_points[0].x() + dt, value(m_points[0].x() + dt));
			else return lerp(t, m_points[N].x() - dt, value(m_points[N].x() - dt), m_points[N].x(), m_points[N].y());
		}
		return 0;
		}
		break;
	case REPEAT:
		{
			if (t < m_points[0].x()) while (t < m_points[0].x()) t += Dt;
			else while (t > m_points[N].x()) t -= Dt;
			return value(t);
		}
		break;
	case REPEAT_OFFSET:
		{
			int n = 0;
			if (t < m_points[0].x()) while (t < m_points[0].x()) { t += Dt; --n; }
			else while (t > m_points[N].x()) { t -= Dt; ++n; }
			double off = n*(m_points[N].y() - m_points[0].y());
			return value(t) + off;
		}
		break;
	}

	return 0;
}

//-----------------------------------------------------------------------------
// This function finds the index of the first load point 
// for which the time is greater than t.
// It returns -1 if t is larger than the last time value
//

int FEPointFunction::FindPoint(double t, double& tval, int startIndex)
{
	switch (m_ext)
	{
	case REPEAT:
	case REPEAT_OFFSET:
	{
		double toff = 0.0;
		while (1)
		{
			double ti = 0;
			for (int i = 0; i<Points(); ++i)
			{
				ti = m_points[i].x() + toff;
				if (ti > t) { tval = ti; return i; }
			}
			toff = ti;
		}
	}
	break;
	default:
		if (startIndex < 0) startIndex = 0;
		if (startIndex >= Points()) return -1;
		for (int i = startIndex; i<Points(); ++i)
		{
			double ti = m_points[i].x();
			if (ti > t) { tval = ti; return i; }
		}
	}
	return -1;
}

//-----------------------------------------------------------------------------

bool FEPointFunction::HasPoint(double t) const
{
	const double tmax = m_points[Points() - 1].x();
	const double eps = 1e-7 * tmax;

	for (int i = 0; i<Points(); ++i) if (fabs(m_points[i].x() - t) < eps) return true;

	return false;
}

//-----------------------------------------------------------------------------

void FEPointFunction::Serialize(DumpStream& ar)
{
	FEFunction1D::Serialize(ar);
	if (ar.IsShallow()) return;

	int n;
	if (ar.IsSaving())
	{
		n = (int)m_fnc; ar << n;
		n = (int)m_ext; ar << n;
		ar << m_points;
	}
	else
	{
		ar >> n; m_fnc = (INTFUNC)n;
		ar >> n; m_ext = (EXTMODE)n;
		ar >> m_points;
	}
}

//-----------------------------------------------------------------------------
double FEPointFunction::derive(double time) const
{
	int N = (int)m_points.size();
	if (N <= 1) return 0;

	double Dt = m_points[N - 1].x() - m_points[0].x();
	double dt = Dt*0.001;
	double t0 = time - dt;
	double t1 = time + dt;

	double v1 = value(t1);
	double v0 = value(t0);

	double D = (v1 - v0) / (2 * dt);

	return D;
}

//-----------------------------------------------------------------------------
FEFunction1D* FEPointFunction::copy()
{
	FEPointFunction* f = new FEPointFunction(GetFEModel());

	f->m_fnc = m_fnc;
	f->m_ext = m_ext;

	Clear();
	for (int i=0; i<Points(); ++i)
	{
		LOADPOINT lp = LoadPoint(i);
		f->Add(lp.time, lp.value);
	}

	return f;
}
