#include "stdafx.h"
#include "FECoreFactory.h"
#include "FECoreBase.h"
#include <assert.h>

//-----------------------------------------------------------------------------
//! constructor
FECoreFactory::FECoreFactory(SUPER_CLASS_ID scid, const char* sztype, int nspec) : m_scid(scid)
{ 
	m_sztype = sztype; 
	m_module = 0; 
	m_spec = nspec;
	m_alloc_id = -1;
}

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

	// store the factory that created the class
	pclass->SetFactoryClass(this);

	// build the class descriptor
	if (pclass->BuildClass() == false) return 0;

	// return the pointer
	return pclass;
}
