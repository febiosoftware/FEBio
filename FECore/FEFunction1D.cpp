#include "stdafx.h"
#include "FEFunction1D.h"
#include "DumpStream.h"

FEFunction1D::FEFunction1D(FEModel* fem) : FECoreBase(fem, FEFUNCTION1D_ID)
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
