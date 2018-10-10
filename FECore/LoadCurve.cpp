#include "stdafx.h"
#include "LoadCurve.h"
#include "DumpStream.h"
#include "FECoreKernel.h"
#include "FEFunction1D.h"

FELoadCurve::FELoadCurve()
{
	m_value = 0;
	m_fnc = nullptr;
}

FELoadCurve::FELoadCurve(FEFunction1D* fnc)
{
	m_value = 0;
	m_fnc = fnc;
}

FELoadCurve::FELoadCurve(const FELoadCurve& lc)
{
	m_value = lc.m_value;
	m_fnc = (lc.m_fnc ? lc.m_fnc->copy() : nullptr);
}

void FELoadCurve::operator = (const FELoadCurve& lc)
{
	m_value = lc.m_value;
}

FELoadCurve::~FELoadCurve()
{
	
}

void FELoadCurve::SetFunction(FEFunction1D* f)
{
	if (m_fnc) delete m_fnc;
	m_fnc = f;
}

void FELoadCurve::Serialize(DumpStream& ar)
{
	if (ar.IsSaving())
	{
		ar << m_value;
		ar << m_fnc->GetTypeStr();
	}
	else
	{
		char szlc[256] = { 0 };
		ar >> m_value;
		ar >> szlc;
		m_fnc = fecore_new<FEFunction1D>(szlc, &ar.GetFEModel());
	}
	m_fnc->Serialize(ar);
}

bool FELoadCurve::CopyFrom(FELoadCurve* lc)
{
	m_value = lc->m_value;
	if (m_fnc) delete m_fnc;
	if (lc->m_fnc) m_fnc = lc->m_fnc->copy();

	return true;
}
