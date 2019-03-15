#pragma once
#include <FECore/FEMeshAdaptor.h>

class FEErosionAdaptor : public FEMeshAdaptor
{
public:
	FEErosionAdaptor(FEModel* fem);

	bool Apply(int iteration) override;

private:
	double	m_maxStress;
	int		m_maxElems;
	int		m_maxIters;
	int		m_metric;

	DECLARE_FECORE_CLASS()
};

class FEMaxStressCriterion : public FEMeshAdaptorCriterion
{
public:
	FEMaxStressCriterion(FEModel* fem);

	bool Check(FEElement& el, double& elemVal) override;

private:
	double	m_maxStress;

	DECLARE_FECORE_CLASS()
};
