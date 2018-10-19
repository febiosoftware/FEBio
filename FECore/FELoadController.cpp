#include "stdafx.h"
#include "FELoadController.h"

FELoadController::FELoadController(FEModel* fem) : FECoreBase(fem, FELOADCONTROLLER_ID)
{
}

void FELoadController::Evaluate(double time)
{
	m_value = GetValue(time);
}

void FELoadController::Serialize(DumpStream& ar)
{
	FECoreBase::Serialize(ar);

	if (ar.IsSaving())
	{
		ar << m_value;
	}
	else
	{
		ar >> m_value;
	}
}
