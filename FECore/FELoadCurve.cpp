/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, Columbia University, and others.

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
#include "FELoadCurve.h"
#include "DumpStream.h"
#include "FECoreKernel.h"
#include "FEFunction1D.h"

BEGIN_FECORE_CLASS(FELoadCurve, FELoadController)
	ADD_PARAMETER(m_fnc.m_points, "point");
	ADD_PARAMETER(m_fnc.m_fnc, "interpolate", FE_PARAM_ATTRIBUTE, "step\0linear\0smooth\0");
	ADD_PARAMETER(m_fnc.m_ext, "extend"     , FE_PARAM_ATTRIBUTE, "constant\0extrapolate\0repeat\0repeat offset\0");
END_FECORE_CLASS();

FELoadCurve::FELoadCurve(FEModel* fem) : FELoadController(fem), m_fnc(fem)
{
}

FELoadCurve::FELoadCurve(const FELoadCurve& lc) : FELoadController(lc), m_fnc(lc.GetFEModel())
{
	m_fnc = lc.m_fnc;
}

void FELoadCurve::operator = (const FELoadCurve& lc)
{
	m_fnc = lc.m_fnc;
}

FELoadCurve::~FELoadCurve()
{
	
}

void FELoadCurve::Serialize(DumpStream& ar)
{
	FELoadController::Serialize(ar);
	m_fnc.Serialize(ar);
}

//! evaluates the loadcurve at time
double FELoadCurve::GetValue(double time)
{
	return m_fnc.value(time);
}

bool FELoadCurve::CopyFrom(FELoadCurve* lc)
{
	m_fnc = lc->m_fnc;
	return true;
}

void FELoadCurve::Add(double time, double value)
{
	m_fnc.Add(time, value);
}

void FELoadCurve::Clear()
{
	m_fnc.Clear();
}

void FELoadCurve::SetInterpolation(FEPointFunction::INTFUNC f)
{
	m_fnc.SetInterpolation(f);
}

void FELoadCurve::SetExtendMode(FEPointFunction::EXTMODE f)
{
	m_fnc.SetExtendMode(f);
}
