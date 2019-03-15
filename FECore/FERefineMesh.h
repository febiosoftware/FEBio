#pragma once
#include "fecore_api.h"
#include "FEMeshAdaptor.h"
class FEModel;

class FEMeshTopo;
class FEFixedBC;
class FEPrescribedDOF;

class FECORE_API FERefineMesh : public FEMeshAdaptor
{
public:
	FERefineMesh(FEModel* fem);

protected:
	bool BuildMeshTopo(FEModel& fem);

	void UpdateFixedBC(FEFixedBC& bc);
	void UpdatePrescribedBC(FEPrescribedDOF& bc);

protected:
	FEMeshTopo*	m_topo;
	int			m_N0;
	int			m_NC;
	int			m_NN;
};
