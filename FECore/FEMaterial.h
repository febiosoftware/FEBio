// FEMaterial.h: interface for the FEMaterial class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FEMATERIAL_H__07F3E572_45B6_444E_A3ED_33FE9D18E82D__INCLUDED_)
#define AFX_FEMATERIAL_H__07F3E572_45B6_444E_A3ED_33FE9D18E82D__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "tens4d.h"
#include "LoadCurve.h"		//---> can we delete this?
#include "FECoordSysMap.h"
#include "FEMaterialFactory.h"
#include "FEParameterList.h"
#include "FEMaterialPoint.h"
#include "DumpFile.h"
#include <string.h>

#define INRANGE(x, a, b) ((x)>=(a) && (x)<=(b))
#define IN_RIGHT_OPEN_RANGE(x, a, b) ((x)>=(a) && (x)<(b))

//-----------------------------------------------------------------------------
//! exception to throw during the material initialization phase

class MaterialError
{
public:
	MaterialError(const char* sz, ...);

	const char* Error() { return m_szerr; }

protected:
	char	m_szerr[512];
};

//-----------------------------------------------------------------------------
//! exception to throw during material initialization phase
class MaterialRangeError
{
public:
	// szvar = name of variable
	// vmin  = inf value
	// vmax  = sup value
	// bl    = inf is allowed
	// br    = sup is allowed
	MaterialRangeError(const char* szvar, double vmin, double vmax, bool bl, bool br) : m_szvar(szvar), m_vmin(vmin), m_vmax(vmax), m_bl(bl), m_br(br) {}

public:
	const char*	m_szvar;
	double	m_vmin, m_vmax;
	bool	m_bl, m_br;
};

//-----------------------------------------------------------------------------
//! Abstract base class for material types

//! From this class all other material classes are derived.

class FEMaterial : public FEParamContainer
{
public:
	FEMaterial();
	virtual ~FEMaterial();

	//! set/get material name
	void SetName(const char* sz) { strcpy(m_szname, sz); }
	const char* GetName() { return m_szname; }

	//! returns a pointer to a new material point object
	virtual FEMaterialPoint* CreateMaterialPointData() { return 0; };

	//! performs initialization and parameter checking
	virtual void Init(){}

	int GetID() { return m_nID; }
	void SetID(int nid) { m_nID = nid; }

	//! Serialize material data to archive
	virtual void Serialize(DumpFile& ar);

protected:
	char	m_szname[128];	//!< name of material
	int		m_nID;			//!< material ID
};

//-----------------------------------------------------------------------------
//! Base class for solid-materials.
//! These materials need to define the stress and tangent functions.
//!
class FESolidMaterial : public FEMaterial
{
public:
	//! calculate stress at material point
	virtual mat3ds Stress(FEMaterialPoint& pt) = 0;

	//! calculate tangent stiffness at material point
	virtual tens4ds Tangent(FEMaterialPoint& pt) = 0;

	//! calculate derivative of stress w.r.t. solute concentration at material point
	mat3ds Tangent_Concentration(FEMaterialPoint& pt, const int isol);
	
	//! return the material density
	virtual double Density() = 0;

	//! return the molar mass
	virtual double MolarMass() = 0;
};

//-----------------------------------------------------------------------------
//! Base class for nested materials. A nested material describes whose material 
//! response depends on the formulation of another user specified material. For
//! instance, the FEViscoElastic is an example of nested materials.
//!
//! NOTE:	The visco-elastic material uses a new format that no longer requires
//!			the m_nBaseMat parameter anymore. In the interest of backward
//!			compatibility I am implementing a hack that allows both formulations.
//!			If the m_nBaseMat is -1, then the new format
//!			is used where the m_pBase is owned by the FENestedMaterial. If
//!			it is not -1, then it uses the old format where m_nBaseMat
//!			points to another material.

class FENestedMaterial : public FESolidMaterial
{
public:
	FENestedMaterial() { m_nBaseMat = -1; m_pBase = 0; }
	virtual ~FENestedMaterial(){}

	//! return solid component's density
	double Density () { return m_pBase->Density(); }

	//! return solid component's molar mass
	double MolarMass () { return m_pBase->MolarMass(); }

	//! serialization
	void Serialize(DumpFile& ar);

public:
	int					m_nBaseMat;	//!< material ID of base material (one-based!)
	FESolidMaterial*	m_pBase;	//!< pointer to base material

	DECLARE_PARAMETER_LIST();
};

#endif // !defined(AFX_FEMATERIAL_H__07F3E572_45B6_444E_A3ED_33FE9D18E82D__INCLUDED_)
