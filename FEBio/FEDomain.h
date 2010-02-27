#pragma once

#include "FEElement.h"
#include "FECore/vector.h"
#include "Archive.h"

class FEMesh;
class FESolidSolver;
class FEHeatSolver;

#define FE_SOLID_DOMAIN			1
#define FE_SHELL_DOMAIN			2
#define FE_SURFACE_DOMAIN		3
#define FE_TRUSS_DOMAIN			4
#define FE_RIGID_SOLID_DOMAIN	5
#define FE_RIGID_SHELL_DOMAIN	6
#define FE_UDGHEX_DOMAIN		7
#define FE_PORO_SOLID_DOMAIN	8
#define FE_HEAT_SOLID_DOMAIN	9

//-----------------------------------------------------------------------------
//! This class describes a physical domain that will be divided into elements
//! of a specific type. All elements in the domain have to have the same type
//!
class FEDomain
{
public:
	FEDomain(FEMesh* pm, int ntype) { m_pMesh = pm; m_ntype = ntype; }
	virtual ~FEDomain() {}

	int Type() { return m_ntype; }

	void SetMesh(FEMesh* pm) { m_pMesh = pm; }

	virtual void create(int n) = 0;

	virtual int Elements() = 0;

	virtual FEDomain* Clone() { assert(false); return 0; }

	virtual void Reset() {}

	virtual void Serialize(FEM& fem, Archive& ar) {}

	virtual bool Initialize(FEM& fem) { return true; }

	// TODO: this is temporary and will be moved to a different class
	virtual void UpdateStresses(FEM& fem) {}

	virtual void InitElements() {}

	virtual void StiffnessMatrix(FESolidSolver* psolver) {}

	virtual void Residual(FESolidSolver* psolver, vector<double>& R) {}

	// TODO: this is not the preferred interface but I've added it for now
	virtual FEElement& ElementRef(int i) = 0;

	FEElement* FindElementFromID(int nid);

	virtual void UnpackElement(FEElement& el, unsigned int nflags = FE_UNPACK_ALL) = 0;

	void SetMatID(int mid)
	{
		for (int i=0; i<Elements(); ++i) ElementRef(i).SetMatID(mid);
	}

protected:
	FEMesh*	m_pMesh;

	int	m_ntype;
};

//-----------------------------------------------------------------------------
//! domain described by Lagrange-type 3D volumetric elements
//!
class FEElasticSolidDomain : public FEDomain
{
public:
	FEElasticSolidDomain(FEMesh* pm) : FEDomain(pm, FE_SOLID_DOMAIN) {}

	void create(int n) { m_Elem.resize(n); }
	int Elements() { return m_Elem.size(); }
	FESolidElement& operator [] (int n) { return m_Elem[n]; }
	
	FESolidElement& Element(int n) { return m_Elem[n]; }

	FEElement& ElementRef(int n) { return m_Elem[n]; }

	FEElasticSolidDomain& operator = (FEElasticSolidDomain& d) { m_Elem = d.m_Elem; m_pMesh = d.m_pMesh; return (*this); }

	FEDomain* Clone()
	{
		FEElasticSolidDomain* pd = new FEElasticSolidDomain(m_pMesh);
		pd->m_Elem = m_Elem; pd->m_pMesh = m_pMesh;
		return pd;
	}

	FEElement* FindElementFromID(int nid);

	bool Initialize(FEM& fem);

	void InitElements();

	void Reset();

	void Serialize(FEM& fem, Archive& ar);

	//! Unpack solid element data
	void UnpackElement(FEElement& el, unsigned int nflag = FE_UNPACK_ALL);

	// update stresses
	void UpdateStresses(FEM& fem);

	//! calculates the global stiffness matrix for this domain
	void StiffnessMatrix(FESolidSolver* psolver);

	//! calculates the solid element stiffness matrix
	void ElementStiffness(FEM& fem, FESolidElement& el, matrix& ke);

	//! calculates the residual
	void Residual(FESolidSolver* psolver, vector<double>& R);

	//! Calculates the internal stress vector for solid elements
	void InternalForces(FESolidElement& el, vector<double>& fe);

protected:
	// --- S T I F F N E S S ---

	//! Dilatational stiffness component for nearly-incompressible materials
	void DilatationalStiffness(FEM& fem, FESolidElement& elem, matrix& ke);

	//! geometrical stiffness (i.e. initial stress)
	void GeometricalStiffness(FESolidElement& el, matrix& ke);

	//! material stiffness component
	void MaterialStiffness(FEM& fem, FESolidElement& el, matrix& ke);

	//! calculates the solid element inertial stiffness matrix
	void ElementInertialStiffness(FEM& fem, FESolidElement& el, matrix& ke);

	// --- R E S I D U A L ---

	//! Calculatess external body forces for solid elements
	void BodyForces(FEM& fem, FESolidElement& elem, vector<double>& fe);

	// ---

protected:
	vector<FESolidElement>	m_Elem;
};

//-----------------------------------------------------------------------------
//! domain class for poro-elastic 3D solid elements
class FEPoroSolidDomain : public FEElasticSolidDomain
{
public:
	//! constructor
	FEPoroSolidDomain(FEMesh* pm) : FEElasticSolidDomain(pm) { m_ntype = FE_PORO_SOLID_DOMAIN; }

	FEDomain* Clone()
	{
		FEPoroSolidDomain* pd = new FEPoroSolidDomain(m_pMesh);
		pd->m_Elem = m_Elem; pd->m_pMesh = m_pMesh;
		return pd;
	}

	//! calculates the global stiffness matrix for this domain
	void StiffnessMatrix(FESolidSolver* psolver);

	//! calculates the residual
	void Residual(FESolidSolver* psolver, vector<double>& R);

	// update stresses
	void UpdateStresses(FEM& fem);

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
	FERigidSolidDomain(FEMesh* pm) : FEElasticSolidDomain(pm) { m_ntype = FE_RIGID_SOLID_DOMAIN; }

	FEDomain* Clone()
	{
		FERigidSolidDomain* pd = new FERigidSolidDomain(m_pMesh);
		pd->m_Elem = m_Elem; pd->m_pMesh = m_pMesh;
		return pd;
	}

	//! calculates the global stiffness matrix for this domain
	void StiffnessMatrix(FESolidSolver* psolver);

	//! calculates the residual
	void Residual(FESolidSolver* psolver, vector<double>& R);

	// update stresses
	void UpdateStresses(FEM& fem);
};

//-----------------------------------------------------------------------------
//! domain class for uniform-deformation-gradient hex elements (UDG)
class FEUDGHexDomain : public FEElasticSolidDomain
{
public:
	//! constructor
	FEUDGHexDomain(FEMesh* pm) : FEElasticSolidDomain(pm) { m_ntype = FE_UDGHEX_DOMAIN; }

	FEDomain* Clone()
	{
		FEUDGHexDomain* pd = new FEUDGHexDomain(m_pMesh);
		pd->m_Elem = m_Elem; pd->m_pMesh = m_pMesh;
		return pd;
	}

	//! calculates the residual
	void Residual(FESolidSolver* psolver, vector<double>& R);

	//! calculates the global stiffness matrix for this domain
	void StiffnessMatrix(FESolidSolver* psolver);

	// update stresses
	void UpdateStresses(FEM& fem);

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

	//! dilatational stiffness for UDG hex elements
	void UDGDilatationalStiffness(FEM& fem, FESolidElement& el, matrix& ke);

	//! geometrical stiffness for UDG hex elements
	void UDGGeometricalStiffness(FEM& fem, FESolidElement& el, matrix& ke);

	//! material stiffness for UDG hex elements
	void UDGMaterialStiffness(FEM& fem, FESolidElement& el, matrix& ke);

protected:
	void AvgCartDerivs(FESolidElement& el, double GX[8], double GY[8], double GZ[8], int state = 0);
	void AvgDefGrad(FESolidElement& el, mat3d& F, double GX[8], double GY[8], double GZ[8]);
	double HexVolume(FESolidElement& el, int state = 0);
};

//-----------------------------------------------------------------------------
//! Domain described by 3D shell elements
class FEElasticShellDomain : public FEDomain
{
public:
	FEElasticShellDomain(FEMesh* pm) : FEDomain(pm, FE_SHELL_DOMAIN) {}

	void create(int n) { m_Elem.resize(n); }
	int Elements() { return m_Elem.size(); }
	FEShellElement& operator [] (int n) { return m_Elem[n]; }

	FEShellElement& Element(int n) { return m_Elem[n]; }

	FEElement& ElementRef(int n) { return m_Elem[n]; }

	FEElasticShellDomain& operator = (FEElasticShellDomain& d) { m_Elem = d.m_Elem; m_pMesh = d.m_pMesh; return (*this); }

	FEDomain* Clone()
	{
		FEElasticShellDomain* pd = new FEElasticShellDomain(m_pMesh);
		pd->m_Elem = m_Elem; pd->m_pMesh = m_pMesh;
		return pd;
	}

	void Reset();

	void InitElements();

	void Serialize(FEM& fem, Archive& ar);

	bool Initialize(FEM& fem);

	//! Unpack shell element data
	void UnpackElement(FEElement& el, unsigned int nflag = FE_UNPACK_ALL);

	// update stresses
	void UpdateStresses(FEM& fem);

	// --- S T I F F N E S S --- 

	//! calculates the global stiffness matrix for this domain
	void StiffnessMatrix(FESolidSolver* psolver);

	//! calculates the shell element stiffness matrix
	void ElementStiffness(FEM& fem, FEShellElement& el, matrix& ke);

	//! Dilatational stiffness component for nearly-incompressible materials
	void DilatationalStiffness(FEM& fem, FEShellElement& elem, matrix& ke);

	// --- R E S I D U A L ---

	//! calculates the residual
	void Residual(FESolidSolver* psolver, vector<double>& R);

	//! Calculates the internal stress vector for shell elements
	void InternalForces(FEShellElement& el, vector<double>& fe);

	//! Calculate extenral body forces for shell elements
	void BodyForces(FEM& fem, FEShellElement& el, vector<double>& fe);

protected:
	vector<FEShellElement>	m_Elem;
};

//-----------------------------------------------------------------------------
//! domain class for 3D rigid shell elements
//!
class FERigidShellDomain : public FEElasticShellDomain
{
public:
	//! constructor
	FERigidShellDomain(FEMesh* pm) : FEElasticShellDomain(pm) { m_ntype = FE_RIGID_SHELL_DOMAIN; }

	FEDomain* Clone()
	{
		FERigidShellDomain* pd = new FERigidShellDomain(m_pMesh);
		pd->m_Elem = m_Elem; pd->m_pMesh = m_pMesh;
		return pd;
	}

	//! calculates the global stiffness matrix for this domain
	void StiffnessMatrix(FESolidSolver* psolver);

	//! calculates the residual
	void Residual(FESolidSolver* psolver, vector<double>& R);

	// update stresses
	void UpdateStresses(FEM& fem);
};

//-----------------------------------------------------------------------------
//! Domain described by 3D truss elements
class FEElasticTrussDomain : public FEDomain
{
public:
	FEElasticTrussDomain(FEMesh* pm) : FEDomain(pm, FE_TRUSS_DOMAIN) {}

	void create(int n) { m_Elem.resize(n); }
	int Elements() { return m_Elem.size(); }
	FETrussElement& operator [] (int n) { return m_Elem[n]; }

	FETrussElement& Element(int i) { return m_Elem[i]; }

	FEElement& ElementRef(int n) { return m_Elem[n]; }

	FEDomain* Clone()
	{
		FEElasticTrussDomain* pd = new FEElasticTrussDomain(m_pMesh);
		pd->m_Elem = m_Elem; pd->m_pMesh = m_pMesh;
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

protected:
	vector<FETrussElement>	m_Elem;
};

//-----------------------------------------------------------------------------
//! domain class for 3D heat elements
class FEHeatSolidDomain : public FEDomain
{
public:
	FEHeatSolidDomain(FEMesh* pm) : FEDomain(pm, FE_HEAT_SOLID_DOMAIN) {}

	void create(int n) { m_Elem.resize(n); }
	int Elements() { return m_Elem.size(); }
	FEDomain* Clone()
	{
		FEHeatSolidDomain* pd = new FEHeatSolidDomain(m_pMesh);
		pd->m_Elem = m_Elem; pd->m_pMesh = m_pMesh;
		return pd;
	}
	FEElement& ElementRef(int n) { return m_Elem[n]; }

	//! Unpack solid element data
	void UnpackElement(FEElement& el, unsigned int nflag = FE_UNPACK_ALL);

	void StiffnessMatrix(FEHeatSolver* psolver);

	FESolidElement& Element(int i) { return m_Elem[i]; }

protected:
	//! calculate the conductive element stiffness matrix
	void ConductionStiffness(FEM& fem, FESolidElement& el, matrix& ke);

	//! calculate the capacitance element stiffness matrix
	void CapacitanceStiffness(FEM& fem, FESolidElement& el, matrix& ke);

protected:
	vector<FESolidElement>	m_Elem;
};
