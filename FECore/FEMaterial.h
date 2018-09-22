// FEMaterial.h: interface for the FEMaterial class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FEMATERIAL_H__07F3E572_45B6_444E_A3ED_33FE9D18E82D__INCLUDED_)
#define AFX_FEMATERIAL_H__07F3E572_45B6_444E_A3ED_33FE9D18E82D__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "tens4d.h"
#include "FECoreBase.h"
#include "FEMaterialPoint.h"
#include "FECoordSysMap.h"
#include "DumpStream.h"
#include "FECoreKernel.h"
#include <string.h>
#include <stddef.h>

#define INRANGE(x, a, b) ((x)>=(a) && (x)<=(b))
#define IN_RIGHT_OPEN_RANGE(x, a, b) ((x)>=(a) && (x)<(b))

//-----------------------------------------------------------------------------
// forward declaration of some classes
class FEModel;
class FEElement;
class FEDomain;

//-----------------------------------------------------------------------------
//! helper functions for reporting material errors

bool FECORE_API MaterialError(const char* sz, ...);

//-----------------------------------------------------------------------------
// Forward declaration of the FEElasticMaterial class. 
// TODO: The only reason I had to do this is to define the FEMaterial::GetElasticMaterial.
// However, this is only a temporary construct so make sure to delete this forward declaration
// when no longer needed.
class FEElasticMaterial;

//-----------------------------------------------------------------------------
//! Abstract base class for material types

//! From this class all other material classes are derived.

class FECORE_API FEMaterial : public FECoreBase
{
public:
	FEMaterial(FEModel* pfem);
	virtual ~FEMaterial();

	//! returns a pointer to a new material point object
	virtual FEMaterialPoint* CreateMaterialPointData() { return 0; };

	//! performs initialization
	bool Init();

	//! Serialize material data to archive
	void Serialize(DumpStream& ar);

	//! Return elastic material \todo I need to move this function up the hierarchy once I redesign the material library
	virtual FEElasticMaterial* GetElasticMaterial() { return 0; }

public:
	//! Set the local coordinate system map
	void SetCoordinateSystemMap(FECoordSysMap* pmap);

	//! Get the local coordinate system
	FECoordSysMap* GetCoordinateSystemMap();

	//! Set the local coordinate for a material point
	virtual void SetLocalCoordinateSystem(FEElement& el, int n, FEMaterialPoint& mp);

	//! return the model this material belongs to
	FEModel* GetFEModel();

public:
	//! Assign a domain to this material
	void AddDomain(FEDomain* dom);


private:
	FECoordSysMap*	m_pmap;			//!< local material coordinate system
	FEModel*		m_pfem;			//!< pointer to model this material belongs to
};

#endif // !defined(AFX_FEMATERIAL_H__07F3E572_45B6_444E_A3ED_33FE9D18E82D__INCLUDED_)
