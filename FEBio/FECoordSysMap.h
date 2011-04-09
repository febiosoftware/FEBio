// FECoordSysMap.h: interface for the FECoordSysMap class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FECOORDSYSMAP_H__5BEAB9FF_6AAE_4CCE_876C_2A2866A8165C__INCLUDED_)
#define AFX_FECOORDSYSMAP_H__5BEAB9FF_6AAE_4CCE_876C_2A2866A8165C__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "mat3d.h"
#include "FEElement.h"
#include "DumpFile.h"

#define FE_MAP_NONE		0
#define FE_MAP_LOCAL	1
#define FE_MAP_SPHERE	2
#define FE_MAP_VECTOR	3

//-----------------------------------------------------------------------------
//! The FECoordSysMap class is used to create local coordinate systems.

class FECoordSysMap  
{
public:
	FECoordSysMap(int ntype) { m_ntype = ntype; }
	virtual ~FECoordSysMap() {}

	//! return the local coordinate system at an element's gauss point
	virtual mat3d LocalElementCoord(FEElement& el, int n) = 0;

	virtual void Serialize(DumpFile& ar) = 0;

public:
	int	m_ntype;
};

//-----------------------------------------------------------------------------

class FELocalMap : public FECoordSysMap
{
public:
	FELocalMap();

	void SetLocalNodes(int n1, int n2, int n3);

	mat3d LocalElementCoord(FEElement& el, int n);

	virtual void Serialize(DumpFile& ar);

protected:
	int	m_n[3];	// local element nodes
};

//-----------------------------------------------------------------------------

class FESphericalMap : public FECoordSysMap
{
public:
	FESphericalMap() : FECoordSysMap(FE_MAP_SPHERE) {}

	void SetSphereCenter(vec3d c) { m_c = c; }

	mat3d LocalElementCoord(FEElement& el, int n);

	virtual void Serialize(DumpFile& ar);

protected:
	vec3d	m_c;	// center of map
};

//-----------------------------------------------------------------------------

class FEVectorMap : public FECoordSysMap
{
public:
	FEVectorMap() : FECoordSysMap(FE_MAP_VECTOR) {}

	void SetVectors(vec3d a, vec3d d) { m_a = a; m_d = d; }

	mat3d LocalElementCoord(FEElement& el, int n);

	virtual void Serialize(DumpFile& ar);

public:
	vec3d	m_a, m_d;
};

#endif // !defined(AFX_FECOORDSYSMAP_H__5BEAB9FF_6AAE_4CCE_876C_2A2866A8165C__INCLUDED_)
