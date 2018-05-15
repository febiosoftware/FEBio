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

//-----------------------------------------------------------------------------
//! helper functions for reporting material errors

bool FECOREDLL_EXPORT MaterialError(const char* sz, ...);

//-----------------------------------------------------------------------------
// Forward declaration of the FEElasticMaterial class. 
// TODO: The only reason I had to do this is to define the FEMaterial::GetElasticMaterial.
// However, this is only a temporary construct so make sure to delete this forward declaration
// when no longer needed.
class FEElasticMaterial;

//-----------------------------------------------------------------------------
//! Abstract base class for material types

//! From this class all other material classes are derived.

class FECOREDLL_EXPORT FEMaterial : public FECoreBase
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

	//! Set the local coordinate for a material point
	virtual void SetLocalCoordinateSystem(FEElement& el, int n, FEMaterialPoint& mp);

	//! return the model this material belongs to
	FEModel* GetFEModel();

private:
	int		m_nRB;			//!< rigid body ID (TODO: I hope to remove this sometime)

private:
	FECoordSysMap*	m_pmap;			//!< local material coordinate system
	FEModel*		m_pfem;			//!< pointer to model this material belongs to
};

#endif // !defined(AFX_FEMATERIAL_H__07F3E572_45B6_444E_A3ED_33FE9D18E82D__INCLUDED_)
