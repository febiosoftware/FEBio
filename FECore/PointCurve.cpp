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
#include "PointCurve.h"
#include "BSpline.h"
#include <assert.h>
#include <algorithm>

#ifndef min
#define min(a,b) ((a)<(b)?(a):(b))
#endif

class PointCurve::Imp
{
public:
	int		fnc;	//!< interpolation function
	int		ext;	//!< extend mode
	std::vector<vec2d>	points;
	BSpline* spline;    //!< B-spline
};

//-----------------------------------------------------------------------------
PointCurve::PointCurve() : im(new PointCurve::Imp)
{
	im->fnc = LINEAR;
	im->ext = CONSTANT;
	im->spline = nullptr;
}

//-----------------------------------------------------------------------------
PointCurve::PointCurve(const PointCurve& pc) : im(new PointCurve::Imp)
{
	im->fnc = pc.im->fnc;
	im->ext = pc.im->ext;
	im->points = pc.im->points;
	im->spline = nullptr;
	Update();
}

//-----------------------------------------------------------------------------
void PointCurve::operator = (const PointCurve& pc)
{
	im->fnc = pc.im->fnc;
	im->ext = pc.im->ext;
	im->points = pc.im->points;
    if (im->spline) delete im->spline;
	Update();
}

//-----------------------------------------------------------------------------
PointCurve::~PointCurve()
{
	if (im->spline) delete im->spline;
}

//-----------------------------------------------------------------------------
//! adds a point to the point curve
int PointCurve::Add(double x, double y)
{
	// find the place to insert the data point
	int n = 0;
	int nsize = Points();
	while ((n < nsize) && (im->points[n].x() < x)) ++n;

	// insert loadpoint
	im->points.insert(im->points.begin() + n, vec2d(x, y));

	return n;
}

//-----------------------------------------------------------------------------
int PointCurve::Add(const vec2d& p)
{
	return Add(p.x(), p.y());
}

//-----------------------------------------------------------------------------
//! Clears the loadcurve data
void PointCurve::Clear()
{
	im->points.clear();
	if (im->spline) delete im->spline;
	im->spline = nullptr;
}

//-----------------------------------------------------------------------------
//! return nr of points
int PointCurve::Points() const
{
	return (int) im->points.size();
}

//-----------------------------------------------------------------------------
// Sets the time and data value of point i
// This function assumes that the load curve data has already been created
//
void PointCurve::SetPoint(int i, double x, double y)
{
	vec2d& pt = im->points[i];
	pt.x() = x;
	pt.y() = y;
}

//-----------------------------------------------------------------------------
void PointCurve::SetPoint(int i, const vec2d& p)
{
	im->points[i] = p;
}

//-----------------------------------------------------------------------------
void PointCurve::SetPoints(const std::vector<vec2d>& points)
{
	im->points = points;
}

//-----------------------------------------------------------------------------
//! return all points
std::vector<vec2d> PointCurve::GetPoints() const
{
	return im->points;
}

//-----------------------------------------------------------------------------
//! Set the type of interpolation
void PointCurve::SetInterpolator(int fnc) 
{ 
	im->fnc = fnc; 
}

//-----------------------------------------------------------------------------
//! return current interpolator
int PointCurve::GetInterpolator() const
{
	return im->fnc;
}

//-----------------------------------------------------------------------------
//! Set the extend mode
void PointCurve::SetExtendMode(int mode) 
{ 
	im->ext = mode; 
}

//-----------------------------------------------------------------------------
//! Get the extend mode
int PointCurve::GetExtendMode() const
{
	return im->ext;
}

//-----------------------------------------------------------------------------
//! get a point
vec2d PointCurve::Point(int i) const
{
	return im->points[i];
}

//-----------------------------------------------------------------------------
void PointCurve::Delete(int n)
{
	if ((n >= 0) && (n < Points()) && (Points() > 2))
	{
		im->points.erase(im->points.begin() + n);
	}
}

//-----------------------------------------------------------------------------
void PointCurve::Delete(const std::vector<int>& indexList)
{
	std::vector<int> tmp;
	int N = (int)indexList.size();
	for (int i = 0; i < N; ++i)
	{
		int n = indexList[i];
		if ((n >= 0) && (n < Points())) tmp.push_back(n);
	}

	std::sort(tmp.begin(), tmp.end());

	for (int i = 0; i < N; ++i)
	{
		int n = tmp[i];
		im->points.erase(im->points.begin() + n);
		for (int j = i + 1; j < N; ++j) tmp[j]--;
	}
}

//-----------------------------------------------------------------------------
void PointCurve::Scale(double s)
{
	for (int i = 0; i < Points(); ++i)
	{
		im->points[i].y() *= s;
	}
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
	return f0 + (f1 - f0) * (t - t0) / (t1 - t0);
}

inline double qerp(double t, double t0, double f0, double t1, double f1, double t2, double f2)
{
	double q0 = ((t2 - t) * (t1 - t)) / ((t2 - t0) * (t1 - t0));
	double q1 = ((t2 - t) * (t - t0)) / ((t2 - t1) * (t1 - t0));
	double q2 = ((t - t1) * (t - t0)) / ((t2 - t1) * (t2 - t0));

	return f0 * q0 + f1 * q1 + f2 * q2;
}

double PointCurve::value(double time) const
{
	std::vector<vec2d>& points = im->points;
	int nsize = Points();
	if (nsize == 0) return 0;
	if (nsize == 1) return points[0].y();

	int N = nsize - 1;

	if (time == points[0].x()) return points[0].y();
	if (time == points[N].x()) return points[N].y();

	double tmax = points[N].x();
	double tmin = points[0].x();
	if ((im->fnc > SMOOTH) && (im->fnc < SMOOTH_STEP))
	{
		if (time > tmax) return im->spline->eval(tmax);
		else if (time < tmin) return im->spline->eval(tmin);
		else return im->spline->eval(time);
	}
    

	if (time < tmin) return ExtendValue(time);
	if (time > tmax) return ExtendValue(time);

	if (im->fnc == LINEAR)
	{
		int n = 0;
		while (points[n].x() <= time) ++n;

		double t0 = points[n - 1].x();
		double t1 = points[n].x();

		double f0 = points[n - 1].y();
		double f1 = points[n].y();

		return lerp(time, t0, f0, t1, f1);
	}
	else if (im->fnc == STEP)
	{
		int n = 0;
		while (points[n].x() <= time) ++n;

		return points[n].y();
	}
	else if (im->fnc == SMOOTH_STEP)
	{
		int n = 0;
		while (points[n].x() <= time) ++n;

		double t0 = points[n - 1].x();
		double t1 = points[n].x();

		double f0 = points[n - 1].y();
		double f1 = points[n].y();

		double w = (time - t0) / (t1 - t0);
		double w2 = w * w;
		double w3 = w * w2;

		return f0 + (f1 - f0)*w3*(10.0 - 15.0*w + 6.0*w2);
	}
	else if (im->fnc == SMOOTH)
	{
		if (nsize == 2)
		{
			double t0 = points[0].x();
			double t1 = points[1].x();

			double f0 = points[0].y();
			double f1 = points[1].y();

			return lerp(time, t0, f0, t1, f1);
		}
		else if (nsize == 3)
		{
			double t0 = points[0].x();
			double t1 = points[1].x();
			double t2 = points[2].x();

			double f0 = points[0].y();
			double f1 = points[1].y();
			double f2 = points[2].y();

			return qerp(time, t0, f0, t1, f1, t2, f2);
		}
		else
		{
			int n = 0;
			while (points[n].x() <= time) ++n;

			if (n == 1)
			{
				double t0 = points[0].x();
				double t1 = points[1].x();
				double t2 = points[2].x();

				double f0 = points[0].y();
				double f1 = points[1].y();
				double f2 = points[2].y();

				return qerp(time, t0, f0, t1, f1, t2, f2);
			}
			else if (n == nsize - 1)
			{
				double t0 = points[n - 2].x();
				double t1 = points[n - 1].x();
				double t2 = points[n].x();

				double f0 = points[n - 2].y();
				double f1 = points[n - 1].y();
				double f2 = points[n].y();

				return qerp(time, t0, f0, t1, f1, t2, f2);
			}
			else
			{
				double t0 = points[n - 2].x();
				double t1 = points[n - 1].x();
				double t2 = points[n].x();
				double t3 = points[n + 1].x();

				double f0 = points[n - 2].y();
				double f1 = points[n - 1].y();
				double f2 = points[n].y();
				double f3 = points[n + 1].y();

				double q1 = qerp(time, t0, f0, t1, f1, t2, f2);
				double q2 = qerp(time, t1, f1, t2, f2, t3, f3);

				return lerp(time, t1, q1, t2, q2);
			}
		}
	}
	return 0;
}

//-----------------------------------------------------------------------------
//! This function determines the value of the point curve outside of its domain
//!
double PointCurve::ExtendValue(double t) const
{
	int nsize = Points();
	int N = nsize - 1;

	std::vector<vec2d>& points = im->points;

	if (nsize == 0) return 0;
	if (nsize == 1) return points[0].y();

	double Dt = (points[N].x() - points[0].x());
	double dt = 0.001 * Dt;
	if (dt == 0) return points[0].y();

	switch (im->ext)
	{
	case CONSTANT:
		if (t < points[0].x()) return points[0].y();
		if (t > points[N].x()) return points[N].y();
		break;
	case EXTRAPOLATE:
		switch (im->fnc)
		{
		case STEP:
		case SMOOTH_STEP:
		{
			if (t < points[0].x()) return points[0].y();
			if (t > points[N].x()) return points[N].y();
		}
		break;
		case LINEAR:
		{
			if (t < points[0].x()) return lerp(t, points[0].x(), points[0].y(), points[1].x(), points[1].y());
			else return lerp(t, points[N - 1].x(), points[N - 1].y(), points[N].x(), points[N].y());
		}
		break;
		case SMOOTH:
		{
			if (t < points[0].x()) return lerp(t, points[0].x(), points[0].y(), points[0].x() + dt, value(points[0].x() + dt));
			else return lerp(t, points[N].x() - dt, value(points[N].x() - dt), points[N].x(), points[N].y());
		}
		return 0;
		}
		break;
	case REPEAT:
	{
		if (t < points[0].x()) while (t < points[0].x()) t += Dt;
		else while (t > points[N].x()) t -= Dt;
		return value(t);
	}
	break;
	case REPEAT_OFFSET:
	{
		int n = 0;
		if (t < points[0].x()) while (t < points[0].x()) { t += Dt; --n; }
		else while (t > points[N].x()) { t -= Dt; ++n; }
		double off = n * (points[N].y() - points[0].y());
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

int PointCurve::FindPoint(double t, double& tval, int startIndex)
{
	switch (im->ext)
	{
	case REPEAT:
	case REPEAT_OFFSET:
	{
		double toff = 0.0;
		while (1)
		{
			double ti = 0;
			for (int i = 0; i < Points(); ++i)
			{
				ti = im->points[i].x() + toff;
				if (ti > t) { tval = ti; return i; }
			}
			toff = ti;
		}
	}
	break;
	default:
		if (startIndex < 0) startIndex = 0;
		if (startIndex >= Points()) return -1;
		for (int i = startIndex; i < Points(); ++i)
		{
			double ti = im->points[i].x();
			if (ti > t) { tval = ti; return i; }
		}
	}
	return -1;
}

//-----------------------------------------------------------------------------

bool PointCurve::HasPoint(double t) const
{
	const double tmax = im->points[Points() - 1].x();
	const double eps = 1e-7 * tmax;

	for (int i = 0; i < Points(); ++i) if (fabs(im->points[i].x() - t) < eps) return true;

	return false;
}

//-----------------------------------------------------------------------------
double PointCurve::derive(double time) const
{
	int N = (int)im->points.size();
	if (N <= 1) return 0;
	double tmax = im->points[N - 1].x();
	double tmin = im->points[0].x();

	if ((im->fnc > SMOOTH) && (im->fnc < SMOOTH_STEP))
	{
		if (time > tmax) return im->spline->eval_deriv(tmax);
		else if (time < tmin) return im->spline->eval_deriv(tmin);
		else return im->spline->eval_deriv(time);
	}

	double Dt = im->points[N - 1].x() - im->points[0].x();
	double dt = Dt * 1e-9;
	double D = 0;

	if (time >= tmax) {
		// use backward difference
		double t2 = time - 2 * dt;
		double t1 = time - dt;
		double t0 = time;

		double v2 = value(t2);
		double v1 = value(t1);
		double v0 = value(t0);

		D = (v2 - 4 * v1 + 3 * v0) / (2 * dt);
	}
	else if (time <= tmin) {
		// use forward difference
		double t0 = time;
		double t1 = time + dt;
		double t2 = time + 2 * dt;

		double v0 = value(t0);
		double v1 = value(t1);
		double v2 = value(t2);

		D = (-v2 + 4 * v1 - 3 * v0) / (2 * dt);
	}
	else {
		// use central difference
		double t0 = time - dt;
		double t1 = time + dt;

		double v1 = value(t1);
		double v0 = value(t0);

		D = (v1 - v0) / (2 * dt);
	}

	return D;
}

//-----------------------------------------------------------------------------
double PointCurve::deriv2(double time) const
{
	int N = (int)im->points.size();
	if (N <= 1) return 0;
	double tmax = im->points[N - 1].x();
	double tmin = im->points[0].x();

	if ((im->fnc > SMOOTH) && (im->fnc < SMOOTH_STEP))
	{
		if (time > tmax) return im->spline->eval_deriv2(tmax);
		else if (time < tmin) return im->spline->eval_deriv2(tmin);
		else return im->spline->eval_deriv2(time);
	}

	double Dt = im->points[N - 1].x() - im->points[0].x();
	double dt = Dt * 1e-3;
	double D = 0;

	if (time >= tmax) {
		// use backward difference
		double t2 = time - 2 * dt;
		double t1 = time - dt;
		double t0 = time;

		double v2 = value(t2);
		double v1 = value(t1);
		double v0 = value(t0);

		D = (v2 - 2 * v1 + v0) / (dt * dt);
	}
	else if (time <= tmin) {
		// use forward difference
		double t0 = time;
		double t1 = time + dt;
		double t2 = time + 2 * dt;

		double v0 = value(t0);
		double v1 = value(t1);
		double v2 = value(t2);

		D = (v2 - 2 * v1 + v0) / (dt * dt);
	}
	else {
		// use central difference
		double t0 = time - dt;
		double t1 = time;
		double t2 = time + dt;

		double v0 = value(t0);
		double v1 = value(t1);
		double v2 = value(t2);

		D = (v2 - 2 * v1 + v0) / (dt * dt);
	}

	return D;
}

double PointCurve::integrate(double a, double b) const
{
	if (a == b) return 0;

	// Swap a and b if a is greater than b
	// negate the answer at the end.
	int neg = 1;
	if (a > b)
	{
		double temp = a;
		a = b;
		b = temp;
		neg = -1;
	}

	double integral = 0.0;

	// if both points are outside the bounds of the load curve 
	// just do a single trapezoid
	// TODO: add cases for  repeat and repeat offset curves
	std::vector<vec2d>& points = im->points;
	if (a > points[Points() - 1].x() || b < points[0].x())
	{
		integral = (b - a) * (value(a) + value(b)) / 2;
	}
	else
	{
		// Find index of first point larger than a
		int start = -1;

		for (int index = 0; index < Points(); index++)
		{
			if (points[index].x() > a)
			{
				start = index;
				break;
			}
		}

		// Do trapezoid rule from a to next point
		integral += (points[start].x() - value(a)) * ((points[start].y() + value(a)) / 2);


		// Loop over points between a and b and do trapezoid rule for each interval
		int index;
		for (index = start; index < Points() - 1; index++)
		{
			// Stop before overshooting b
			if (points[index + 1].x() >= b) break;
			integral += (points[index + 1].x() - points[index].x()) * ((points[index + 1].y() + points[index].y()) * 0.5);
		}

		//Do trapezoid rule from most recent point to b
		integral += (value(b) - points[index].x()) * ((value(b) + points[index].y()) * 0.5);
	}

	return integral * neg;
}

bool PointCurve::Update()
{
	bool bvalid = true;

	if ((im->fnc > SMOOTH) && (im->fnc < SMOOTH_STEP))
	{
		const int N = Points();

		// initialize B-spline
        if (im->spline) delete im->spline;
		im->spline = new BSpline();
		switch (im->fnc) {
		case CSPLINE:
		{
			int korder = min(N, 4);
			if (!im->spline->init_interpolation(korder, im->points)) bvalid = false;
		}
		break;
		case CPOINTS:
		{
			int korder = min(N, 4);
			if (!im->spline->init(korder, im->points)) bvalid = false;
		}
		break;
		case APPROX:
		{
			int korder = min(N / 2 + 1, 4);
			if (!im->spline->init_approximation(korder, N / 2 + 1, im->points)) bvalid = false;
		}
		break;
		default:
			bvalid = false;
			assert(false);
		}

		if (bvalid == false) { delete im->spline; im->spline = nullptr; }
	}
	return bvalid;
}
