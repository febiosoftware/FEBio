#include "stdafx.h"
#include "FELoadController.h"

REGISTER_SUPER_CLASS(FELoadController, FELOADCONTROLLER_ID);

FELoadController::FELoadController(FEModel* fem) : FECoreBase(fem)
{
}

void FELoadController::Evaluate(double time)
{
	m_value = GetValue(time);
}

void FELoadController::Serialize(DumpStream& ar)
{
	FECoreBase::Serialize(ar);
	ar & m_value;
}
