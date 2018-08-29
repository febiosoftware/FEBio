#include "stdafx.h"
#include "FECoreFactory.h"
#include <assert.h>

//-----------------------------------------------------------------------------
//! constructor
FECoreFactory::FECoreFactory(SUPER_CLASS_ID scid, const char* sztype) : m_scid(scid) { m_sztype = sztype; m_module = 0; }

//-----------------------------------------------------------------------------
//! virtual constructor
FECoreFactory::~FECoreFactory(){}

//-----------------------------------------------------------------------------
//! set the module ID
void FECoreFactory::SetModuleID(unsigned int mid)
{
	m_module = mid;
}

//-----------------------------------------------------------------------------
void* FECoreFactory::CreateInstance(FEModel* pfem)
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
