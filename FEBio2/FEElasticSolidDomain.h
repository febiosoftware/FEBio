#pragma once
#include "FECore/FESolidDomain.h"

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
	void StiffnessMatrix(FESolidSolver* psolver);

	//! calculates the solid element stiffness matrix
	virtual void ElementStiffness(FEM& fem, int iel, matrix& ke);

	//! calculates the residual
	void Residual(FESolidSolver* psolver, vector<double>& R);

	//! Calculates the internal stress vector for solid elements
	void InternalForces(FESolidElement& el, vector<double>& fe);

protected:
	// --- S T I F F N E S S ---

	//! geometrical stiffness (i.e. initial stress)
	void GeometricalStiffness(FESolidElement& el, matrix& ke);

	//! material stiffness component
	void MaterialStiffness(FEM& fem, FESolidElement& el, matrix& ke);

	//! calculates the solid element inertial stiffness matrix
	void ElementInertialStiffness(FEM& fem, FESolidElement& el, matrix& ke);

	//! calculates the stiffness matrix due to body forces 
	void BodyForceStiffness(FEM& fem, FESolidElement& el, matrix& ke);

	// --- R E S I D U A L ---

	//! Calculatess external body forces for solid elements
	void BodyForces(FEM& fem, FESolidElement& elem, vector<double>& fe);

	// ---
};

//-----------------------------------------------------------------------------
//! The following domain implements the finite element formulation for a three-field
//! volume element. 
class FE3FieldElasticSolidDomain : public FEElasticSolidDomain
{
protected:
	struct ELEM_DATA
	{
		double	eJ;		// average element jacobian
		double	ep;		// average pressure
		double	Lk;		// Lagrangian multiplier
	};

public:
	//! constructor
	FE3FieldElasticSolidDomain(FEMesh* pm, FEMaterial* pmat) : FEElasticSolidDomain(pm, pmat) { m_ntype = FE_3F_SOLID_DOMAIN; }

	//! TODO: do I really use this?
	FE3FieldElasticSolidDomain& operator = (FE3FieldElasticSolidDomain& d) { m_Elem = d.m_Elem; m_pMesh = d.m_pMesh; return (*this); }

	//! create a clone of this class
	FEDomain* Clone()
	{
		FE3FieldElasticSolidDomain* pd = new FE3FieldElasticSolidDomain(m_pMesh, m_pMat);
		pd->m_Elem = m_Elem; pd->m_pMesh = m_pMesh; pd->m_Node = m_Node;
		pd->m_Data = m_Data;
		return pd;
	}

	//! initialize class
	bool Initialize(FEModel& fem);

	// update stresses
	void UpdateStresses(FEModel& fem);

	//! calculates the solid element stiffness matrix
	void ElementStiffness(FEM& fem, int iel, matrix& ke);

	//! augmentation
	bool Augment();

	//! serialize data to archive
	void Serialize(DumpFile& ar);

protected:
	//! Dilatational stiffness component for nearly-incompressible materials
	void DilatationalStiffness(FEM& fem, int iel, matrix& ke);

	//! material stiffness component
	void MaterialStiffness(FEM& fem, int iel, matrix& ke);

	//! geometrical stiffness (i.e. initial stress)
	void GeometricalStiffness(int iel, matrix& ke);

protected:
	vector<ELEM_DATA>	m_Data;
};

//-----------------------------------------------------------------------------
//! Domain class for poro-elastic 3D solid elements
//! Note that this class inherits from FEElasticSolidDomain since the poro domain
//! also needs to calculate elastic stiffness contributions.
//!
class FEPoroSolidDomain : public FEElasticSolidDomain
{
public:
	//! constructor
	FEPoroSolidDomain(FEMesh* pm, FEMaterial* pmat) : FEElasticSolidDomain(pm, pmat) { m_ntype = FE_PORO_SOLID_DOMAIN; }

	FEDomain* Clone()
	{
		FEPoroSolidDomain* pd = new FEPoroSolidDomain(m_pMesh, m_pMat);
		pd->m_Elem = m_Elem; pd->m_pMesh = m_pMesh; pd->m_Node = m_Node;
		return pd;
	}

	//! calculates the global stiffness matrix for this domain
	void StiffnessMatrix(FESolidSolver* psolver);

	//! calculates the residual
	void Residual(FESolidSolver* psolver, vector<double>& R);

	// update stresses
	void UpdateStresses(FEModel& fem);

protected:
	//! Calculates the internal fluid forces
	bool InternalFluidWork(FEM& fem, FESolidElement& elem, vector<double>& fe);

	//! calculates the element biphasic stiffness matrix
	bool ElementPoroStiffness(FEM& fem, FESolidElement& el, matrix& ke);

	//! calculates the solid element stiffness matrix
	void SolidElementStiffness(FEM& fem, FESolidElement& el, matrix& ke);

	//! material stiffness component
	void PoroMaterialStiffness(FEM& fem, FESolidElement& el, matrix& ke);
};

//-----------------------------------------------------------------------------
//! domain class for 3D rigid elements
//!
class FERigidSolidDomain : public FEElasticSolidDomain
{
public:
	//! constructor
	FERigidSolidDomain(FEMesh* pm, FEMaterial* pmat) : FEElasticSolidDomain(pm, pmat) { m_ntype = FE_RIGID_SOLID_DOMAIN; }

	FEDomain* Clone()
	{
		FERigidSolidDomain* pd = new FERigidSolidDomain(m_pMesh, m_pMat);
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

//-----------------------------------------------------------------------------
//! domain class for uniform-deformation-gradient hex elements (UDG)
class FEUDGHexDomain : public FEElasticSolidDomain
{
public:
	//! constructor
	FEUDGHexDomain(FEMesh* pm, FEMaterial* pmat) : FEElasticSolidDomain(pm, pmat) { m_ntype = FE_UDGHEX_DOMAIN; }

	FEDomain* Clone()
	{
		FEUDGHexDomain* pd = new FEUDGHexDomain(m_pMesh, m_pMat);
		pd->m_Elem = m_Elem; pd->m_pMesh = m_pMesh; pd->m_Node = m_Node;
		return pd;
	}

	//! calculates the residual
	void Residual(FESolidSolver* psolver, vector<double>& R);

	//! calculates the global stiffness matrix for this domain
	void StiffnessMatrix(FESolidSolver* psolver);

	// update stresses
	void UpdateStresses(FEModel& fem);

protected:
	//! Calculates the internal stress vector for enhanced strain hex elements
	void UDGInternalForces(FEM& fem, FESolidElement& el, vector<double>& fe);

	//! calculates hourglass forces for the UDG element
	void UDGHourglassForces(FEM& fem, FESolidElement& el, vector<double>& fe);

protected:
	//! calculate element stiffness
	void UDGElementStiffness(FEM& fem, FESolidElement& el, matrix& ke);

	//! hourglass stiffness for UDG hex elements
	void UDGHourglassStiffness(FEM& fem, FESolidElement& el, matrix& ke);

	//! geometrical stiffness for UDG hex elements
	void UDGGeometricalStiffness(FEM& fem, FESolidElement& el, matrix& ke);

	//! material stiffness for UDG hex elements
	void UDGMaterialStiffness(FEM& fem, FESolidElement& el, matrix& ke);

protected:
	void AvgCartDerivs(FESolidElement& el, double GX[8], double GY[8], double GZ[8], int state = 0);
	void AvgDefGrad(FESolidElement& el, mat3d& F, double GX[8], double GY[8], double GZ[8]);
	double HexVolume(FESolidElement& el, int state = 0);
};
