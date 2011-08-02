#pragma once
#include "FECore/FEShellDomain.h"

//-----------------------------------------------------------------------------
// forward declaration of FEM class
class FEM;

//-----------------------------------------------------------------------------
//! Domain described by 3D shell elements
class FEElasticShellDomain : public FEShellDomain
{
public:
	FEElasticShellDomain(FEMesh* pm, FEMaterial* pmat) : FEShellDomain(FE_SHELL_DOMAIN, pm, pmat) {}

	//! TODO: do I really need this?
	FEElasticShellDomain& operator = (FEElasticShellDomain& d) { m_Elem = d.m_Elem; m_pMesh = d.m_pMesh; return (*this); }

	FEDomain* Clone()
	{
		FEElasticShellDomain* pd = new FEElasticShellDomain(m_pMesh, m_pMat);
		pd->m_Elem = m_Elem; pd->m_pMesh = m_pMesh; pd->m_Node = m_Node;
		return pd;
	}

	void Reset();

	void InitElements();

	bool Initialize(FEModel& fem);

	void Serialize(DumpFile& ar);

	//! Unpack shell element data
	void UnpackLM(FEElement& el, vector<int>& lm);

	// update stresses
	void UpdateStresses(FEModel& fem);

	// --- S T I F F N E S S --- 

	//! calculates the global stiffness matrix for this domain
	void StiffnessMatrix(FESolidSolver* psolver);

	//! calculates the shell element stiffness matrix
	void ElementStiffness(FEM& fem, int iel, matrix& ke);

	// --- R E S I D U A L ---

	//! calculates the residual
	void Residual(FESolidSolver* psolver, vector<double>& R);

	//! Calculates the internal stress vector for shell elements
	void InternalForces(FEShellElement& el, vector<double>& fe);

	//! Calculate extenral body forces for shell elements
	void BodyForces(FEM& fem, FEShellElement& el, vector<double>& fe);
};

//-----------------------------------------------------------------------------
//! domain class for 3D rigid shell elements
//!
class FERigidShellDomain : public FEElasticShellDomain
{
public:
	//! constructor
	FERigidShellDomain(FEMesh* pm, FEMaterial* pmat) : FEElasticShellDomain(pm, pmat) { m_ntype = FE_RIGID_SHELL_DOMAIN; }

	FEDomain* Clone()
	{
		FERigidShellDomain* pd = new FERigidShellDomain(m_pMesh, m_pMat);
		pd->m_Elem = m_Elem; pd->m_pMesh = m_pMesh; pd->m_Node = m_Node;
		return pd;
	}

	//! calculates the global stiffness matrix for this domain
	void StiffnessMatrix(FESolidSolver* psolver);

	//! calculates the residual
	void Residual(FESolidSolver* psolver, vector<double>& R);

	// update stresses
	void UpdateStresses(FEModel& fem);
};
