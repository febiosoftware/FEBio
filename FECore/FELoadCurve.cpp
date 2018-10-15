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
	m_fnc = nullptr;
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
