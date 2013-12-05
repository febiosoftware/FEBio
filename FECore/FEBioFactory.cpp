#include "stdafx.h"
#include "FEBioFactory.h"
#include <assert.h>

//-----------------------------------------------------------------------------
//! constructor
FEBioFactory::FEBioFactory(SUPER_CLASS_ID scid, const char* sztype) : m_scid(scid) { m_sztype = sztype; }

//-----------------------------------------------------------------------------
//! virtual constructor
FEBioFactory::~FEBioFactory(){}

//-----------------------------------------------------------------------------
void* FEBioFactory::CreateInstance(FEModel* pfem)
{
	// create a new instance of this class
	FECoreBase* pclass = static_cast<FECoreBase*>(Create(pfem)); assert(pclass);
	if (pclass == 0) return 0;

	// make sure the super-class ID matches
	// TODO: I need to delete pclass
	if (pclass->GetSuperClassID() != GetSuperClassID()) return 0;

	// set the type string of this class
	pclass->SetTypeStr(GetTypeStr());

	// return the pointer
	return pclass;
}
