#pragma once
#include "FEMeshAdaptor.h"

class FEErosionAdaptor : public FEMeshAdaptor
{
public:
	FEErosionAdaptor(FEModel* fem);

	bool Apply(int iteration) override;

private:
	int		m_maxIters;

	FEMeshAdaptorCriterion*	m_criterion;

	DECLARE_FECORE_CLASS()
};
