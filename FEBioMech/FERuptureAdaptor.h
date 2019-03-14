#pragma once
#include <FECore/FEMeshAdaptor.h>

class FERuptureAdaptor : public FEMeshAdaptor
{
public:
	FERuptureAdaptor(FEModel* fem);

	bool Apply() override;

private:
	double	m_maxStress;
	int		m_maxElems;
	int		m_maxIters;
	int		m_metric;

	DECLARE_FECORE_CLASS()
};
