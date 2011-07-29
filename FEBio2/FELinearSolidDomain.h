#pragma once

#include "FECore/FESolidDomain.h"

//-----------------------------------------------------------------------------
// Forward declarations
class FEM;
class FELinearSolidSolver;

//-----------------------------------------------------------------------------
//! Class describing a linear elastic solid domain
class FELinearSolidDomain : public FESolidDomain
{
public:
	FELinearSolidDomain(FEMesh* pm, FEMaterial* pmat) : FESolidDomain(FE_LINEAR_SOLID_DOMAIN, pm, pmat) {}

	FEDomain* Clone()
	{
		FELinearSolidDomain* pd = new FELinearSolidDomain(m_pMesh, m_pMat);
		pd->m_Elem = m_Elem; pd->m_pMesh = m_pMesh; pd->m_Node = m_Node;
		return pd;
	}

	//! Unpack solid element data
	void UnpackLM(FEElement& el, vector<int>& lm);

	void StiffnessMatrix(FELinearSolidSolver* psolver);

	void UpdateStresses(FEModel& fem);

protected:
	void ElementStiffness(FEM& fem, FESolidElement& el, matrix& ke);
};
