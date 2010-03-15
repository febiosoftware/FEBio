#pragma once

#include "FEElement.h"
#include "FECore/vector.h"
#include "Archive.h"

class FEMesh;
class FESolidSolver;
class FEHeatSolver;
class FEMaterial;

#define FE_SOLID_DOMAIN			1
#define FE_SHELL_DOMAIN			2
#define FE_SURFACE_DOMAIN		3
#define FE_TRUSS_DOMAIN			4
#define FE_RIGID_SOLID_DOMAIN	5
#define FE_RIGID_SHELL_DOMAIN	6
#define FE_UDGHEX_DOMAIN		7
#define FE_UT4_DOMAIN			8
#define FE_PORO_SOLID_DOMAIN	9
#define FE_HEAT_SOLID_DOMAIN	10

//-----------------------------------------------------------------------------
//! This class describes a physical domain that will be divided into elements
//! of a specific type. All elements in the domain have to have the same type
//!
class FEDomain
{
public:
	FEDomain(int ntype, FEMesh* pm, FEMaterial* pmat) { m_pMesh = pm; m_ntype = ntype; m_pMat = pmat; }
	virtual ~FEDomain() {}

	int Type() { return m_ntype; }

	void SetMesh(FEMesh* pm) { m_pMesh = pm; }
	FEMesh* GetMesh() { return m_pMesh; }

	void SetMaterial(FEMaterial* pmat) { m_pMat = pmat; }
	FEMaterial* GetMaterial() { return m_pMat; }

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

	//!< Initialize material point data for the elements
	void InitMaterialPointData();

	// TODO: this is not the preferred interface but I've added it for now
	virtual FEElement& ElementRef(int i) = 0;

	FEElement* FindElementFromID(int nid);

	virtual void UnpackElement(FEElement& el, unsigned int nflags = FE_UNPACK_ALL) = 0;

	void SetMatID(int mid)
	{
		for (int i=0; i<Elements(); ++i) ElementRef(i).SetMatID(mid);
	}

protected:
	FEMesh*		m_pMesh;	//!< the mesh that this domain is a part of
	FEMaterial*	m_pMat;		//!< the material for this domain

	int	m_ntype;
};

//-----------------------------------------------------------------------------
//! abstract base class for 3D volumetric elements
class FESolidDomain : public FEDomain
{
public:
	//! constructor
	FESolidDomain(int ntype, FEMesh* pm, FEMaterial* pmat) : FEDomain(ntype, pm, pmat) {}

	//! create storage for elements
	void create(int nsize) { m_Elem.resize(nsize); }

	//! return nr of elements
	int Elements() { return m_Elem.size(); }

	//! element access
	FESolidElement& Element(int n) { return m_Elem[n]; }
	FESolidElement& ElementRef(int n) { return m_Elem[n]; }

protected:
	vector<FESolidElement>	m_Elem;	//!< array of elements
};

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
		pd->m_Elem = m_Elem; pd->m_pMesh = m_pMesh;
		return pd;
	}

	//! initialize class
	bool Initialize(FEM& fem);

	//! Init elements
	void InitElements();

	//! reset element data
	void Reset();

	//! serialize data to archive
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
	FERigidSolidDomain(FEMesh* pm, FEMaterial* pmat) : FEElasticSolidDomain(pm, pmat) { m_ntype = FE_RIGID_SOLID_DOMAIN; }

	FEDomain* Clone()
	{
		FERigidSolidDomain* pd = new FERigidSolidDomain(m_pMesh, m_pMat);
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
	FEUDGHexDomain(FEMesh* pm, FEMaterial* pmat) : FEElasticSolidDomain(pm, pmat) { m_ntype = FE_UDGHEX_DOMAIN; }

	FEDomain* Clone()
	{
		FEUDGHexDomain* pd = new FEUDGHexDomain(m_pMesh, m_pMat);
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
//! Abstract base class for shell elements
class FEShellDomain : public FEDomain
{
public:
	//! constructor
	FEShellDomain(int ntype, FEMesh* pm, FEMaterial* pmat) : FEDomain(ntype, pm, pmat) {}

	//! create storage for elements
	void create(int nsize) { m_Elem.resize(nsize); }

	//! return nr of elements
	int Elements() { return m_Elem.size(); }

	//! element access
	FEShellElement& Element(int n) { return m_Elem[n]; }
	FEElement& ElementRef(int n) { return m_Elem[n]; }

protected:
	vector<FEShellElement>	m_Elem;	//!< array of elements
};

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
//! Abstract base class for truss elements
class FETrussDomain : public FEDomain
{
public:
	FETrussDomain(int ntype, FEMesh* pm, FEMaterial* pmat) : FEDomain(ntype, pm, pmat){}

	void create(int n) { m_Elem.resize(n); }
	int Elements() { return m_Elem.size(); }

	FETrussElement& Element(int i) { return m_Elem[i]; }

	FEElement& ElementRef(int n) { return m_Elem[n]; }

protected:
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
};

//-----------------------------------------------------------------------------
//! domain class for 3D heat elements
class FEHeatSolidDomain : public FESolidDomain
{
public:
	FEHeatSolidDomain(FEMesh* pm, FEMaterial* pmat) : FESolidDomain(FE_HEAT_SOLID_DOMAIN, pm, pmat) {}

	FEDomain* Clone()
	{
		FEHeatSolidDomain* pd = new FEHeatSolidDomain(m_pMesh, m_pMat);
		pd->m_Elem = m_Elem; pd->m_pMesh = m_pMesh;
		return pd;
	}

	//! Unpack solid element data
	void UnpackElement(FEElement& el, unsigned int nflag = FE_UNPACK_ALL);

	void HeatStiffnessMatrix(FEHeatSolver* psolver);

protected:
	//! calculate the conductive element stiffness matrix
	void ConductionStiffness(FEM& fem, FESolidElement& el, matrix& ke);

	//! calculate the capacitance element stiffness matrix
	void CapacitanceStiffness(FEM& fem, FESolidElement& el, matrix& ke);
};
