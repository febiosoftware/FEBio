#pragma once
#include "FERefineMesh.h"

class FEHexRefine2D : public FERefineMesh
{
public:
	FEHexRefine2D(FEModel* fem);

	bool Apply(int iteration) override;

protected:
	bool RefineMesh(FEModel& fem);
	bool BuildSplitLists(FEModel& fem);
	void UpdateNewNodes(FEModel& fem);
	void FindHangingNodes(FEModel& fem);
	void BuildNewDomains(FEModel& fem);

private:
	int		m_maxelem;			// max nr of elements
	int		m_elemRefine;		// max nr of elements to refine per step
	int		m_maxiter;
	vector<int>	m_elemList;

	int	m_splitElems;
	int m_splitFaces;
	int m_splitEdges;
	int	m_hangingNodes;

	FEMeshAdaptorCriterion*	m_criterion;

	DECLARE_FECORE_CLASS();
};
