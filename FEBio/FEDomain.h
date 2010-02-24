#pragma once

#include "FEElement.h"
#include "FECore/vector.h"
#include "Archive.h"

class FEMesh;

//-----------------------------------------------------------------------------
//! This class describes a physical domain that will be divided into elements
//! of a specific type. All elements in the domain have to have the same type
//!
class FEDomain
{
public:
	FEDomain(FEMesh* pm) { m_pMesh = pm; }
	virtual ~FEDomain() {}

	void SetMesh(FEMesh* pm) { m_pMesh = pm; }

	virtual void create(int n) = 0;

	virtual int Elements() = 0;

	virtual void Reset() {}

	virtual void Serialize(FEM& fem, Archive& ar) {}

protected:
	FEMesh*	m_pMesh;
};

//-----------------------------------------------------------------------------
//! domain described by Lagrange-type 3D volumetric elements
//!
class FESolidDomain : public FEDomain
{
public:
	FESolidDomain(FEMesh* pm) : FEDomain(pm) {}

	void create(int n) { m_Elem.resize(n); }
	int Elements() { return m_Elem.size(); }
	FESolidElement& operator [] (int n) { return m_Elem[n]; }
	
	FESolidElement& Element(int n) { return m_Elem[n]; }

	FESolidDomain& operator = (FESolidDomain& d) { m_Elem = d.m_Elem; m_pMesh = d.m_pMesh; return (*this); }

	FEElement* FindElementFromID(int nid);

	bool Init(FEM& fem);

	void Reset();

	void Serialize(FEM& fem, Archive& ar);

	//! Unpack solid element data
	void UnpackElement(FESolidElement& el, unsigned int nflag = FE_UNPACK_ALL);

	// update stresses
	void UpdateStresses(FEM& fem);

	// --- S T I F F N E S S ---

	//! calculates the solid element stiffness matrix
	void ElementStiffness(FEM& fem, FESolidElement& el, matrix& ke);

	//! Dilatational stiffness component for nearly-incompressible materials
	void DilatationalStiffness(FEM& fem, FESolidElement& elem, matrix& ke);

	//! geometrical stiffness (i.e. initial stress)
	void GeometricalStiffness(FESolidElement& el, matrix& ke);

	//! material stiffness component
	void MaterialStiffness(FEM& fem, FESolidElement& el, matrix& ke);

	//! calculates the solid element inertial stiffness matrix
	void ElementInertialStiffness(FEM& fem, FESolidElement& el, matrix& ke);

	//! calculates the element biphasic stiffness matrix
	bool ElementPoroStiffness(FEM& fem, FESolidElement& el, matrix& ke);

	//! hourglass stiffness for UDG hex elements
	void UDGHourglassStiffness(FEM& fem, FESolidElement& el, matrix& ke);

	//! dilatational stiffness for UDG hex elements
	void UDGDilatationalStiffness(FEM& fem, FESolidElement& el, matrix& ke);

	//! geometrical stiffness for UDG hex elements
	void UDGGeometricalStiffness(FEM& fem, FESolidElement& el, matrix& ke);

	//! material stiffness for UDG hex elements
	void UDGMaterialStiffness(FEM& fem, FESolidElement& el, matrix& ke);

	// --- R E S I D U A L ---

	//! Calculates the internal stress vector for solid elements
	void InternalForces(FESolidElement& el, vector<double>& fe);

	//! Calculates the internal stress vector for enhanced strain hex elements
	void UDGInternalForces(FEM& fem, FESolidElement& el, vector<double>& fe);

	//! calculates hourglass forces for the UDG element
	void UDGHourglassForces(FEM& fem, FESolidElement& el, vector<double>& fe);

	//! Calculatess external body forces for solid elements
	void BodyForces(FEM& fem, FESolidElement& elem, vector<double>& fe);

	//! Calculates the internal fluid forces
	bool InternalFluidWork(FEM& fem, FESolidElement& elem, vector<double>& fe);

	// ---

	void AvgCartDerivs(FESolidElement& el, double GX[8], double GY[8], double GZ[8], int state = 0);
	void AvgDefGrad(FESolidElement& el, mat3d& F, double GX[8], double GY[8], double GZ[8]);
	double HexVolume(FESolidElement& el, int state = 0);

protected:
	vector<FESolidElement>	m_Elem;
};

//-----------------------------------------------------------------------------
//! Domain described by 3D shell elements
class FEShellDomain : public FEDomain
{
public:
	FEShellDomain(FEMesh* pm) : FEDomain(pm) {}

	void create(int n) { m_Elem.resize(n); }
	int Elements() { return m_Elem.size(); }
	FEShellElement& operator [] (int n) { return m_Elem[n]; }

	FEShellElement& Element(int n) { return m_Elem[n]; }

	FEShellDomain& operator = (FEShellDomain& d) { m_Elem = d.m_Elem; m_pMesh = d.m_pMesh; return (*this); }

	FEElement* FindElementFromID(int nid);

	void Reset();

	void Serialize(FEM& fem, Archive& ar);

	bool Init(FEM& fem);

	//! Unpack shell element data
	void UnpackElement(FEShellElement& el, unsigned int nflag = FE_UNPACK_ALL);

	// update stresses
	void UpdateStresses(FEM& fem);

	// --- S T I F F N E S S --- 

	//! calculates the shell element stiffness matrix
	void ElementStiffness(FEM& fem, FEShellElement& el, matrix& ke);

	//! Dilatational stiffness component for nearly-incompressible materials
	void DilatationalStiffness(FEM& fem, FEShellElement& elem, matrix& ke);

	// --- R E S I D U A L ---

	//! Calculates the internal stress vector for shell elements
	void InternalForces(FEShellElement& el, vector<double>& fe);

	//! Calculate extenral body forces for shell elements
	void BodyForces(FEM& fem, FEShellElement& el, vector<double>& fe);

protected:
	vector<FEShellElement>	m_Elem;
};

//-----------------------------------------------------------------------------
//! Domain described by 3D truss elements
class FETrussDomain : public FEDomain
{
public:
	FETrussDomain(FEMesh* pm) : FEDomain(pm) {}

	void create(int n) { m_Elem.resize(n); }
	int Elements() { return m_Elem.size(); }
	FETrussElement& operator [] (int n) { return m_Elem[n]; }

	FETrussElement& Element(int i) { return m_Elem[i]; }

	FETrussDomain& operator = (FETrussDomain& d) { m_Elem = d.m_Elem; m_pMesh = d.m_pMesh; return (*this); }

	FEElement* FindElementFromID(int nid);

	bool Init(FEM& fem) { return true; }

	void Reset();

	//! Unpack truss element data
	void UnpackElement(FETrussElement& el, unsigned int flag = FE_UNPACK_ALL);

	//! calculates the truss element stiffness matrix
	void ElementStiffness(FEM& fem, FETrussElement& el, matrix& ke);

	//! Calculates the internal stress vector for solid elements
	void InternalForces(FETrussElement& el, vector<double>& fe);

	//! update the truss stresses
	void UpdateStresses(FEM& fem);

protected:
	vector<FETrussElement>	m_Elem;
};
