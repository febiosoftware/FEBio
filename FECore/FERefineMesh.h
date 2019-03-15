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

	void DoMeshReset(bool b);

	bool Apply(int iteration) override;

protected:
	bool DoMeshRefinement(FEModel& fem);
	bool DoTetRefinement(FEModel& fem);
	bool DoHexRefinement(FEModel& fem);
	bool BuildMeshTopo(FEModel& fem);

	void UpdateFixedBC(FEFixedBC& bc);
	void UpdatePrescribedBC(FEPrescribedDOF& bc);

private:
	bool		m_doMeshReset;
	int			m_maxelem;

private:
	FEMeshTopo*	m_topo;
	int			m_N0;
	int			m_NC;
	int			m_NN;
	vector<int>	m_tag;

	DECLARE_FECORE_CLASS();
};
