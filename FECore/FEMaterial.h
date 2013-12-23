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
#include "FECoreBase.h"
#include "FEMaterialPoint.h"
#include "FECoordSysMap.h"
#include "DumpFile.h"
#include <string.h>

#define INRANGE(x, a, b) ((x)>=(a) && (x)<=(b))
#define IN_RIGHT_OPEN_RANGE(x, a, b) ((x)>=(a) && (x)<(b))

//-----------------------------------------------------------------------------
// forward declaration of FEModel class
class FEModel;

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
// Forward declaration of the FEElasticMaterial class. 
// TODO: The only reason I had to do this is to define the FEMaterial::GetElasticMaterial.
// However, this is only a temporary construct so make sure to delete this forward declaration
// when no longer needed.
class FEElasticMaterial;

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

class FEMaterial : public FECoreBase
{
public:
	FEMaterial(FEModel* pfem);
	virtual ~FEMaterial();

	//! set material name
	void SetName(const char* sz);

	//! get the material's name
	const char* GetName();

	//! returns a pointer to a new material point object
	virtual FEMaterialPoint* CreateMaterialPointData() { return 0; };

	//! performs initialization and parameter checking
	virtual void Init();

	int GetID() { return m_nID; }
	void SetID(int nid) { m_nID = nid; }

	//! Serialize material data to archive
	virtual void Serialize(DumpFile& ar);

	//! Return elastic material \todo I need to move this function up the hierarchy once I redesign the material library
	virtual FEElasticMaterial* GetElasticMaterial() { return 0; }

public:
	// TODO: Some rigid body stuff is moved to here to avoid RTTI and the definition of rigid materials in FECore, 
	//       as well as simplify some initialization. I hope someday to refactor this a bit.
	//! is this a rigid material
	virtual bool IsRigid() { return false; }

	//! get the ID of the rigid body this material is assigned to (-1 if not)
	int GetRigidBodyID() { return m_nRB; }

	//! Set the rigid body ID this material is assigned to
	void SetRigidBodyID(int rid) { m_nRB = rid; }

	//! return the density
	//! TODO: This was added here because the rigid bodies need it to determine the COM
	virtual double Density() { return 0.0; }

public:
	//! Set the local coordinate system map
	void SetCoordinateSystemMap(FECoordSysMap* pmap);

	//! Get the local coordinate system
	FECoordSysMap* GetCoordinateSystemMap();

	//! return the model this material belongs to
	FEModel* GetFEModel();

	//! Get the parent of this material (zero if none)
	FEMaterial* GetParent() { return m_pParent; }

	//! Set the parent of this material
	void SetParent(FEMaterial* pmat) { m_pParent = pmat; }

public: // interface for getting/setting material properties

	//! get the number of material properties
	virtual int Properties () { return 0; }

	//! get a specific material property
	virtual FEMaterial* GetProperty(int i) { return 0; }

	//! find a material property index ( returns <0 for error)
	virtual int FindPropertyIndex(const char* szname) { return -1; }

	//! set a material property (returns false on error)
	virtual bool SetProperty(int i, FEMaterial* pm) { return false; }

private:
	char	m_szname[128];	//!< name of material
	int		m_nID;			//!< material ID
	int		m_nRB;			//!< rigid body ID (TODO: I hope to remove this sometime)

private:
	FECoordSysMap*	m_pmap;			//!< local material coordinate system
	FEModel*		m_pfem;			//!< pointer to model this material belongs to
	FEMaterial*		m_pParent;		//!< pointer to "parent" material (if any)
};

#endif // !defined(AFX_FEMATERIAL_H__07F3E572_45B6_444E_A3ED_33FE9D18E82D__INCLUDED_)
