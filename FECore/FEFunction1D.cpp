#include "stdafx.h"
#include "FEFunction1D.h"
#include "FEModel.h"
#include "LoadCurve.h"
#include "DumpStream.h"

FEFunction1D::FEFunction1D(FEModel* pfem) : m_fem(*pfem)
{
	m_nlc = -1;
	m_scale = 0.0;
}

void FEFunction1D::SetLoadCurveIndex(int lc, double scale)
{
	m_nlc = lc;
	m_scale = scale;
}

double FEFunction1D::value(double x) const
{
	if (m_nlc < 0) return m_scale;

	FELoadCurve* plc = m_fem.GetLoadCurve(m_nlc);
	if (plc == 0) return m_scale;

	return m_scale*plc->Value(x);
}

double FEFunction1D::derive(double x) const
{
	if (m_nlc < 0) return 0.0;

	FELoadCurve* plc = m_fem.GetLoadCurve(m_nlc);
	if (plc == 0) return 0.0;

	return m_scale*plc->Deriv(x);
}

void FEFunction1D::Serialize(DumpStream& ar)
{
	if (ar.IsSaving())
	{
		ar << m_nlc << m_scale;
	}
	else
	{
		ar >> m_nlc >> m_scale;
	}
}
