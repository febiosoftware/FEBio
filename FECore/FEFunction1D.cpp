#include "stdafx.h"
#include "FEFunction1D.h"
#include "DumpStream.h"

REGISTER_SUPER_CLASS(FEFunction1D, FEFUNCTION1D_ID);

FEFunction1D::FEFunction1D(FEModel* fem) : FECoreBase(fem)
{
}

double FEFunction1D::derive(double x) const
{
	const double eps = 1e-6;
	return (value(x + eps) - value(x))/eps;
}

void FEFunction1D::Serialize(DumpStream& ar)
{
	FECoreBase::Serialize(ar);

	if (ar.IsSaving())
	{
	}
	else
	{
	}
}
