#pragma once
#include <FECore/FEMeshAdaptor.h>


class FEMaxStressCriterion : public FEMeshAdaptorCriterion
{
public:
	FEMaxStressCriterion(FEModel* fem);

	bool Check(FEElement& el, double& elemVal) override;

private:
	double	m_maxStress;
	int		m_metric;

	DECLARE_FECORE_CLASS()
};

