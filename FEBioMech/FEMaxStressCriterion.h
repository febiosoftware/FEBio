#pragma once
#include <FECore/FEMeshAdaptor.h>


class FEMaxStressCriterion : public FEMeshAdaptorCriterion
{
public:
	FEMaxStressCriterion(FEModel* fem);

	bool Check(FEElement& el, double& elemVal) override;

private:
	double	m_maxStress;

	DECLARE_FECORE_CLASS()
};

