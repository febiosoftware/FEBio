#pragma once
#include "FEDomain.h"

//-----------------------------------------------------------------------------
//! domain for discrete elements
class FEDiscreteDomain : public FEDomain
{
public:
	FEDiscreteDomain(int ntype, FEMesh* pm, FEMaterial* pmat) : FEDomain(ntype, pm, pmat) {}

	void create(int n) { m_Elem.resize(n); }
	int Elements() { return (int) m_Elem.size(); }
	FEElement& ElementRef(int n) { return m_Elem[n]; }

	bool Initialize(FEModel& fem);

	int Nodes() { return (int) m_Node.size(); }
	FENode& Node(int i);

protected:
	vector<int>					m_Node;
	vector<FEDiscreteElement>	m_Elem;
};
