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
#include "FEMaterialFactory.h"
#include "FEParameterList.h"
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

class FEMaterial : public FEParamContainer
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

	//! is this a rigid material \todo this is temporary solution to avoid RTTI and the need to define rigid materials in FECore
	virtual bool IsRigid() { return false; }

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

public: // interface for managing attributes

	//! Set the attribute
	virtual bool SetAttribute(const char* szname, const char* szval) { return true; }

private:
	char	m_szname[128];	//!< name of material
	int		m_nID;			//!< material ID

private:
	FECoordSysMap*	m_pmap;			//!< local material coordinate system
	FEModel*		m_pfem;			//!< pointer to model this material belongs to
	FEMaterial*		m_pParent;		//!< pointer to "parent" material (if any)
};

//-----------------------------------------------------------------------------
//! Global solute data
//! This structure uniquely identifies a solute in multiphasic problems
//! \todo Move this to a different file
class FESoluteData 
{
public:
	FESoluteData() { m_nID = -1; m_rhoT = 1; m_M = 1; m_z = 0; m_szname[0] = 0; }

public:
	int		m_nID;			//!< solute ID
	double	m_rhoT;			//!< true solute density
	double	m_M;			//!< solute molecular weight
	int		m_z;			//!< solute charge number
	char	m_szname[128];	//!< solute name
};


//-----------------------------------------------------------------------------
//! Global solid-bound molecule (SBM) data.
//! \todo move this to a different file
class FESBMData 
{
public:
	FESBMData() { m_nID = -1; m_rhoT = 1; m_M = 1; m_z = 0; m_szname[0] = 0; }
	
public:
	int		m_nID;			//!< SBM ID
	double	m_rhoT;			//!< SBM true density
	double	m_M;			//!< SBM molar mass
	int		m_z;			//!< SBM charge number
	char	m_szname[128];	//!< SBM name
};


#endif // !defined(AFX_FEMATERIAL_H__07F3E572_45B6_444E_A3ED_33FE9D18E82D__INCLUDED_)
