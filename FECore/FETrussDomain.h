#pragma once
#include "FEDomain.h"

//-----------------------------------------------------------------------------
//! Abstract base class for truss elements
class FETrussDomain : public FEDomain
{
public:
	FETrussDomain(int ntype, FEMesh* pm, FEMaterial* pmat) : FEDomain(ntype, pm, pmat){}

	void create(int n) { m_Elem.resize(n); }
	int Elements() { return m_Elem.size(); }

	FETrussElement& Element(int i) { return m_Elem[i]; }

	FEElement& ElementRef(int n) { return m_Elem[n]; }

	bool Initialize(FEModel& fem);

	int Nodes() { return (int) m_Node.size(); }
	FENode& Node(int i);

protected:
	vector<int>				m_Node;
	vector<FETrussElement>	m_Elem;
};
