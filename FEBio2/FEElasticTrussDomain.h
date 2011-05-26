#pragma once
#include "FECore/FETrussDomain.h"

//-----------------------------------------------------------------------------
//! Domain described by 3D truss elements
class FEElasticTrussDomain : public FETrussDomain
{
public:
	FEElasticTrussDomain(FEMesh* pm, FEMaterial* pmat) : FETrussDomain(FE_TRUSS_DOMAIN, pm, pmat) {}

	FEDomain* Clone()
	{
		FEElasticTrussDomain* pd = new FEElasticTrussDomain(m_pMesh, m_pMat);
		pd->m_Elem = m_Elem; pd->m_pMesh = m_pMesh; pd->m_Node = m_Node;
		return pd;
	}

	FEElasticTrussDomain& operator = (FEElasticTrussDomain& d) { m_Elem = d.m_Elem; m_pMesh = d.m_pMesh; return (*this); }

	void Reset();

	void InitElements();

	//! Unpack truss element data
	void UnpackLM(FEElement& el, vector<int>& lm);

	//! calculates the global stiffness matrix for this domain
	void StiffnessMatrix(FESolidSolver* psolver);

	//! calculates the truss element stiffness matrix
	void ElementStiffness(FEM& fem, FETrussElement& el, matrix& ke);

	//! calculates the residual
	void Residual(FESolidSolver* psolver, vector<double>& R);

	//! Calculates the internal stress vector for solid elements
	void InternalForces(FETrussElement& el, vector<double>& fe);

	//! update the truss stresses
	void UpdateStresses(FEM& fem);
};
