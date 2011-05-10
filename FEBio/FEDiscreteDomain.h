#pragma once
#include "FEDomain.h"

//-----------------------------------------------------------------------------
//! domain for discrete elements
class FEDiscreteDomain : public FEDomain
{
public:
	FEDiscreteDomain(FEMesh* pm, FEMaterial* pmat) : FEDomain(FE_DISCRETE_DOMAIN, pm, pmat) {}

	FEDomain* Clone()
	{
		FEDiscreteDomain* pd = new FEDiscreteDomain(m_pMesh, m_pMat);
		pd->m_Elem = m_Elem; pd->m_pMesh = m_pMesh; pd->m_Node = m_Node;
		return pd;
	}

	void UnpackElement(FEElement& el, unsigned int nflag = FE_UNPACK_ALL);

	void create(int n) { m_Elem.resize(n); }
	int Elements() { return (int) m_Elem.size(); }
	FEElement& ElementRef(int n) { return m_Elem[n]; }

	void StiffnessMatrix(FESolidSolver* psolver);

	void Residual(FESolidSolver* psolver, vector<double>& R);

	bool Initialize(FEM& fem);

	int Nodes() { return (int) m_Node.size(); }
	FENode& Node(int i);

	void Serialize(DumpFile& ar);

protected:
	vector<int>					m_Node;
	vector<FEDiscreteElement>	m_Elem;
};
