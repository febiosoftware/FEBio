// FECoordSysMap.h: interface for the FECoordSysMap class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FECOORDSYSMAP_H__5BEAB9FF_6AAE_4CCE_876C_2A2866A8165C__INCLUDED_)
#define AFX_FECOORDSYSMAP_H__5BEAB9FF_6AAE_4CCE_876C_2A2866A8165C__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "mat3d.h"
#include "FECoreBase.h"

//-----------------------------------------------------------------------------
// Forward declarations
class FEModel;
class FEMesh;
class FEElement;

//-----------------------------------------------------------------------------
//! The FECoordSysMap class is used to create local coordinate systems.

class FECoordSysMap : public FECoreBase
{
public:
	FECoordSysMap(FEModel* pfem);
	virtual ~FECoordSysMap();

	//! initialization
	bool Init();

	//! return the local coordinate system at an element's gauss point
	virtual mat3d LocalElementCoord(FEElement& el, int n) = 0;

	//! serialization
	virtual void Serialize(DumpStream& ar) = 0;

public:
	FEModel* GetFEModel() { return m_pfem; }

private:
	FEModel*	m_pfem;
};

//-----------------------------------------------------------------------------
//! This class generates a material axes based on the local element node numbering.
class FELocalMap : public FECoordSysMap
{
public:
	FELocalMap(FEModel* pfem);

	bool Init() override;

	void SetLocalNodes(int n1, int n2, int n3);

	mat3d LocalElementCoord(FEElement& el, int n) override;

	virtual void Serialize(DumpStream& ar) override;

public:
	int			m_n[3];	// local element nodes

protected:
	DECLARE_PARAMETER_LIST();
};

//-----------------------------------------------------------------------------
//! This class generates material axes based on a spherical map. 
class FESphericalMap : public FECoordSysMap
{
public:
	FESphericalMap(FEModel* pfem);

	bool Init() override;

	void SetSphereCenter(const vec3d& c) { m_c = c; }

	void SetSphereVector(const vec3d& r) { m_r = r;}

	mat3d LocalElementCoord(FEElement& el, int n) override;

	virtual void Serialize(DumpStream& ar) override;

public:
	vec3d		m_c;	// center of map
	vec3d		m_r;	// vector for parallel transport

protected:
	DECLARE_PARAMETER_LIST();
};


//-----------------------------------------------------------------------------

class FECylindricalMap : public FECoordSysMap
{
public:
	FECylindricalMap(FEModel* pfem);

	bool Init() override;

	void SetCylinderCenter(vec3d c) { m_c = c; }

	void SetCylinderAxis(vec3d a) { m_a = a; m_a.unit(); }

	void SetCylinderRef(vec3d r) { m_r = r; m_r.unit(); }

	mat3d LocalElementCoord(FEElement& el, int n) override;

	virtual void Serialize(DumpStream& ar) override;

public:
	vec3d		m_c;	// center of map
	vec3d		m_a;	// axis
	vec3d		m_r;	// reference direction

protected:
	DECLARE_PARAMETER_LIST();
};

//-----------------------------------------------------------------------------

class FEPolarMap : public FECoordSysMap
{
public:
	FEPolarMap(FEModel* pfem);

	bool Init() override;

	void SetCenter(vec3d c) { m_c = c; }

	void SetAxis(vec3d a) { m_a = a; m_a.unit(); }

	void SetVector0(vec3d r) { m_d0 = r; m_d0.unit(); }
	void SetVector1(vec3d r) { m_d1 = r; m_d1.unit(); }

	void SetRadius0(double r) { m_R0 = r; }
	void SetRadius1(double r) { m_R1 = r; }

	mat3d LocalElementCoord(FEElement& el, int n) override;

	virtual void Serialize(DumpStream& ar) override;

public:
	vec3d		m_c;		// center of map
	vec3d		m_a;		// axis
	vec3d		m_d0, m_d1;	// reference direction
	double		m_R0, m_R1;

protected:
	DECLARE_PARAMETER_LIST();
};

//-----------------------------------------------------------------------------

class FEVectorMap : public FECoordSysMap
{
public:
	FEVectorMap(FEModel* pfem);

	bool Init() override;

	void SetVectors(vec3d a, vec3d d);

	mat3d LocalElementCoord(FEElement& el, int n) override;

	virtual void Serialize(DumpStream& ar) override;

public:
	vec3d	m_a, m_d;

	DECLARE_PARAMETER_LIST();
};

//-----------------------------------------------------------------------------
class FESphericalAngleMap : public FECoordSysMap
{
public:
	FESphericalAngleMap(FEModel* pfem);

	bool Init() override;

	void SetAngles(double theta, double phi);

	mat3d LocalElementCoord(FEElement& el, int n) override;

	virtual void Serialize(DumpStream& ar) override;

public:
	double	m_theta;
	double	m_phi;

	DECLARE_PARAMETER_LIST();
};

#endif // !defined(AFX_FECOORDSYSMAP_H__5BEAB9FF_6AAE_4CCE_876C_2A2866A8165C__INCLUDED_)
