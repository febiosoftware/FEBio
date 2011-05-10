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

	bool Initialize(FEM& fem);

	int Nodes() { return (int) m_Node.size(); }
	FENode& Node(int i);

protected:
	vector<int>				m_Node;
	vector<FETrussElement>	m_Elem;
};

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

	bool Initialize(FEM& fem) { return true; }

	void Reset();

	void InitElements();

	//! Unpack truss element data
	void UnpackElement(FEElement& el, unsigned int flag = FE_UNPACK_ALL);

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
