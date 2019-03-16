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
	bool BuildMeshTopo();

	void UpdateBCs();

private:
	void UpdateFixedBC(FEFixedBC& bc);
	void UpdatePrescribedBC(FEPrescribedDOF& bc);

protected:
	FEMeshTopo*	m_topo;
	int			m_N0;
	int			m_NC;
	int			m_NN;
	vector<int>	m_edgeList;	// list of edge flags to see whether the edge was split
	vector<int>	m_faceList;	// list of face flags to see whether the face was split
};
