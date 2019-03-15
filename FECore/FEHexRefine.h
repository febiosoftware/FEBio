#pragma once
#include "FERefineMesh.h"

class FEHexRefine : public FERefineMesh
{
public:
	FEHexRefine(FEModel* fem);

	bool Apply(int iteration) override;

protected:
	bool DoHexRefinement(FEModel& fem);

	bool RefineMesh(FEModel& fem);
	void UpdateBCs(FEModel& fem);

private:
	int		m_maxelem;
	int		m_maxiter;

	FEMeshAdaptorCriterion*	m_criterion;

	DECLARE_FECORE_CLASS();
};
