#pragma once
#include "FECore/FEDiscreteDomain.h"

//-----------------------------------------------------------------------------
//! domain for discrete elements
class FEDiscreteSpringDomain : public FEDiscreteDomain
{
public:
	FEDiscreteSpringDomain(FEMesh* pm, FEMaterial* pmat) : FEDiscreteDomain(FE_DISCRETE_DOMAIN, pm, pmat) {}

	FEDomain* Clone()
	{
		FEDiscreteSpringDomain* pd = new FEDiscreteSpringDomain(m_pMesh, m_pMat);
		pd->m_Elem = m_Elem; pd->m_pMesh = m_pMesh; pd->m_Node = m_Node;
		return pd;
	}

	void UnpackElement(FEElement& el, unsigned int nflag = FE_UNPACK_ALL);
	void UnpackLM(FEElement& el, vector<int>& lm);

	void StiffnessMatrix(FESolidSolver* psolver);

	void Residual(FESolidSolver* psolver, vector<double>& R);

	void Serialize(DumpFile& ar);
};
