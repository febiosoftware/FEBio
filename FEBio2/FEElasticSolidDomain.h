#pragma once
#include "FECore/FESolidDomain.h"
#include "FECore/FEModel.h"

class FEM;

//-----------------------------------------------------------------------------
//! domain described by Lagrange-type 3D volumetric elements
//!
class FEElasticSolidDomain : public FESolidDomain
{
public:
	//! constructor
	FEElasticSolidDomain(FEMesh* pm, FEMaterial* pmat) : FESolidDomain(FE_SOLID_DOMAIN, pm, pmat) {}

	//! TODO: do I really use this?
	FEElasticSolidDomain& operator = (FEElasticSolidDomain& d) { m_Elem = d.m_Elem; m_pMesh = d.m_pMesh; return (*this); }

	//! create a clone of this class
	FEDomain* Clone()
	{
		FEElasticSolidDomain* pd = new FEElasticSolidDomain(m_pMesh, m_pMat);
		pd->m_Elem = m_Elem; pd->m_pMesh = m_pMesh; pd->m_Node = m_Node;
		return pd;
	}

	//! serialize data to archive
	void Serialize(DumpFile& ar);

	//! initialize class
	bool Initialize(FEModel& fem);

	//! Init elements
	void InitElements();

	//! reset element data
	void Reset();

	//! Unpack solid element data
	void UnpackLM(FEElement& el, vector<int>& lm);

	// update stresses
	void UpdateStresses(FEModel& fem);

	//! calculates the global stiffness matrix for this domain
	void StiffnessMatrix(FENLSolver* psolver);

	//! calculates the solid element stiffness matrix
	virtual void ElementStiffness(FEM& fem, int iel, matrix& ke);

	//! calculates the residual
	void Residual(FENLSolver* psolver, vector<double>& R);

	//! Calculates the internal stress vector for solid elements
	void InternalForces(FESolidElement& el, vector<double>& fe);

protected:
	// --- S T I F F N E S S ---

	//! geometrical stiffness (i.e. initial stress)
	void GeometricalStiffness(FESolidElement& el, matrix& ke);

	//! material stiffness component
	void MaterialStiffness(FEM& fem, FESolidElement& el, matrix& ke);

	//! calculates the solid element inertial stiffness matrix
	void ElementInertialStiffness(FEModel& fem, FESolidElement& el, matrix& ke);

	//! calculates the stiffness matrix due to body forces 
	void BodyForceStiffness(FEModel& fem, FESolidElement& el, matrix& ke);

	// --- R E S I D U A L ---

	//! Calculatess external body forces for solid elements
	void BodyForces(FEModel& fem, FESolidElement& elem, vector<double>& fe);

	// ---
};
