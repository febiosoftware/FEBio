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
#include "FEPointFunction.h"
#include "DumpStream.h"
#include "log.h"
#include "BSpline.h"

#ifndef min
#define min(a,b) ((a)<(b)?(a):(b))
#endif

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEPointFunction, FEFunction1D)
	ADD_PARAMETER(m_int, "interpolate", 0, "linear\0step\0smooth\0cubic spline\0control points\0approximation\0");
	ADD_PARAMETER(m_ext, "extend"     , 0, "constant\0extrapolate\0repeat\0repeat offset\0");
    ADD_PARAMETER(m_bln, "log")->SetFlags(FE_PARAM_HIDDEN);
	ADD_PARAMETER(m_points, "points");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! default constructor
FEPointFunction::FEPointFunction(FEModel* fem) : FEFunction1D(fem)
{
	m_int = PointCurve::LINEAR;
	m_ext = PointCurve::CONSTANT;
    m_bln = false;
}

//-----------------------------------------------------------------------------
FEPointFunction::~FEPointFunction()
{

}

//-----------------------------------------------------------------------------
//! Clears the loadcurve data
bool FEPointFunction::Init()
{
	m_fnc.SetInterpolator(m_int);
	m_fnc.SetExtendMode(m_ext);
	m_fnc.SetPoints(m_points);
	if (m_fnc.Update() == false) return false;

    return FEFunction1D::Init();
}

//-----------------------------------------------------------------------------
//! Clears the loadcurve data
void FEPointFunction::Clear()
{ 
	m_fnc.Clear();
}

//-----------------------------------------------------------------------------
//! return nr of points
int FEPointFunction::Points() const
{ 
	return (int) m_points.size(); 
}

//-----------------------------------------------------------------------------
//! set the points
void FEPointFunction::SetPoints(const std::vector<vec2d>& pts)
{
	m_points = pts;
}

//-----------------------------------------------------------------------------
// Sets the time and data value of point i
// This function assumes that the load curve data has already been created
//
void FEPointFunction::SetPoint(int i, double x, double y)
{
	m_fnc.SetPoint(i, x, y);
}

//-----------------------------------------------------------------------------
//! Set the type of interpolation
void FEPointFunction::SetInterpolation(int fnc) { m_int = fnc; }

//-----------------------------------------------------------------------------
//! Set the extend mode
void FEPointFunction::SetExtendMode(int mode) { m_ext = mode; }

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
	int nsize = m_points.size();
	while ((n<nsize) && (m_points[n].x() < x)) ++n;

	// insert loadpoint
	m_points.insert(m_points.begin() + n, vec2d(x, y));
}

//-----------------------------------------------------------------------------
void FEPointFunction::Scale(double s)
{
	for (int i = 0; i < Points(); ++i)
	{
		m_points[i].y() *= s;
	}
}

//-----------------------------------------------------------------------------
double FEPointFunction::value(double time) const
{
    if (m_bln) time = (time > 0) ? log(time) : m_points[0].x();
	return m_fnc.value(time);
}

//-----------------------------------------------------------------------------
void FEPointFunction::Serialize(DumpStream& ar)
{
	FEFunction1D::Serialize(ar);
	if (ar.IsShallow()) return;

	if (ar.IsSaving())
	{
		ar << m_int << m_ext;
		ar << m_points;
	}
	else
	{
		ar >> m_int >> m_ext;
		ar >> m_points;

		m_fnc.Clear();
		m_fnc.SetInterpolator(m_int);
		m_fnc.SetExtendMode(m_ext);
		m_fnc.SetPoints(m_points);
		m_fnc.Update();
	}
}

//-----------------------------------------------------------------------------
double FEPointFunction::derive(double time) const
{
    if (m_bln) time = (time > 0) ? log(time) : m_points[0].x();
	return m_fnc.derive(time);
}

//-----------------------------------------------------------------------------
double FEPointFunction::deriv2(double time) const
{
    if (m_bln) time = (time > 0) ? log(time) : m_points[0].x();
	return m_fnc.deriv2(time);
}

double FEPointFunction::integrate(double a, double b) const
{
	return m_fnc.integrate(a, b);
}

//-----------------------------------------------------------------------------
FEFunction1D* FEPointFunction::copy()
{
	FEPointFunction* f = new FEPointFunction(GetFEModel());

	f->m_int = m_int;
	f->m_ext = m_ext;
	f->m_points = m_points;
	f->m_fnc = m_fnc;
	return f;
}

//-----------------------------------------------------------------------------
void FEPointFunction::CopyFrom(const FEPointFunction& f)
{
	m_int = f.m_int;
	m_ext = f.m_ext;
	m_points = f.m_points;
	m_fnc = f.m_fnc;
}

//-----------------------------------------------------------------------------
void FEPointFunction::CopyFrom(const PointCurve& f)
{
	m_int = f.GetInterpolator();
	m_ext = f.GetExtendMode();
	m_points = f.GetPoints();
	m_fnc = f;
}
