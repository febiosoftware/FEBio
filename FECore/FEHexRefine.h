#pragma once
#include "FERefineMesh.h"

class FEHexRefine : public FERefineMesh
{
public:
	FEHexRefine(FEModel* fem);

	bool Apply(int iteration) override;

protected:
	bool DoHexRefinement(FEModel& fem);

	void RefineMesh(FEModel& fem);
	void UpdateBCs(FEModel& fem);

private:
	int			m_maxelem;

	DECLARE_FECORE_CLASS();
};
