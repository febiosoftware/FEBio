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
	void BuildSplitLists(FEModel& fem);
	void UpdateNewNodes(FEModel& fem);
	void FindHangingNodes(FEModel& fem);
	void BuildNewDomains(FEModel& fem);
	void UpdateBCs(FEModel& fem);

private:
	int		m_maxelem;
	int		m_maxiter;
	vector<int>	m_elemList;

	int	m_splitElems;
	int m_splitFaces;
	int m_splitEdges;

	FEMeshAdaptorCriterion*	m_criterion;

	DECLARE_FECORE_CLASS();
};
