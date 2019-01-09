#pragma once
#include "fecore_api.h"
class FEModel;

class FEMeshTopo;

class FECORE_API FERefineMesh
{
public:
	FERefineMesh();

	bool Apply(FEModel& fem);

protected:
	bool DoMeshRefinement(FEModel& fem);
	bool DoTetRefinement(FEModel& fem);
	bool DoHexRefinement(FEModel& fem);
	bool BuildMeshTopo(FEModel& fem);

private:
	FEMeshTopo*	m_topo;
};
