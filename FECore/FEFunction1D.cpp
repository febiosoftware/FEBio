#include "stdafx.h"
#include "FEFunction1D.h"
#include "DumpStream.h"

FEFunction1D::FEFunction1D(FEModel* pfem) : FECoreBase(FEFUNCTION1D_ID), m_fem(pfem)
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
