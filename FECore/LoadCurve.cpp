#include "stdafx.h"
#include "LoadCurve.h"
#include "DumpStream.h"

FELoadCurve::FELoadCurve() : FECoreBase(FELOADCURVE_ID)
{
	m_value = 0;
}

FELoadCurve::FELoadCurve(const FELoadCurve& lc) : FECoreBase(FELOADCURVE_ID)
{
	m_value = lc.m_value;
}

void FELoadCurve::operator = (const FELoadCurve& lc)
{
	m_value = lc.m_value;
}

FELoadCurve::~FELoadCurve()
{
	
}

void FELoadCurve::Serialize(DumpStream& ar)
{
	if (ar.IsSaving())
	{
		ar << m_value;
	}
	else
	{
		ar >> m_value;
	}
}
