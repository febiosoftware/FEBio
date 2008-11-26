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

//-----------------------------------------------------------------------------
//! The FECoordSysMap class is used to create local coordinate systems.

class FECoordSysMap  
{
public:
	FECoordSysMap() {}
	virtual ~FECoordSysMap() {}

	//! return the local coordinate system at an element's gauss point
	virtual mat3d LocalElementCoord(FEElement& el, int n) = 0;
};

//-----------------------------------------------------------------------------

class FELocalMap : public FECoordSysMap
{
public:
	FELocalMap();

	void SetLocalNodes(int n1, int n2, int n3);

	mat3d LocalElementCoord(FEElement& el, int n);

protected:
	int	m_n[3];	// local element nodes
};

//-----------------------------------------------------------------------------

class FESphericalMap : public FECoordSysMap
{
public:
	FESphericalMap(){}

	void SetSphereCenter(vec3d c) { m_c = c; }

	mat3d LocalElementCoord(FEElement& el, int n);

protected:
	vec3d	m_c;	// center of map
};

//-----------------------------------------------------------------------------

class FEVectorMap : public FECoordSysMap
{
public:
	FEVectorMap(){}

	void SetVectors(vec3d a, vec3d d) { m_a = a; m_d = d; }

	mat3d LocalElementCoord(FEElement& el, int n);

public:
	vec3d	m_a, m_d;
};

//-----------------------------------------------------------------------------

class FERandom2DMap : public FECoordSysMap
{
public:
	FERandom2DMap(){}

	mat3d LocalElementCoord(FEElement& el, int n);
};


#endif // !defined(AFX_FECOORDSYSMAP_H__5BEAB9FF_6AAE_4CCE_876C_2A2866A8165C__INCLUDED_)
