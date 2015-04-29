#include "stdafx.h"
#include "FEModelComponent.h"
#include <string.h>

//-----------------------------------------------------------------------------
//! The constructor takes two arguments: the SUPER_CLASS_ID which defines the 
//! type of the model component and a pointer to the FEModel object this component
//! belongs to.
FEModelComponent::FEModelComponent(SUPER_CLASS_ID sid, FEModel* pfem) : FECoreBase(sid)
{
	m_pfem = pfem;
	m_bactive = true;
	m_szname = 0;
}

//-----------------------------------------------------------------------------
FEModelComponent::~FEModelComponent()
{
	if (m_szname) delete [] m_szname;
}

//-----------------------------------------------------------------------------
FEModel* FEModelComponent::GetFEModel()
{
	return m_pfem;
}

//-----------------------------------------------------------------------------
const char* FEModelComponent::GetName()
{
	return m_szname;
}

//-----------------------------------------------------------------------------
void FEModelComponent::SetName(const char* sz)
{
	if (sz == 0) return;
	int l = strlen(sz);
	if (m_szname) delete [] m_szname;
	m_szname = new char[l+1];
	if (l > 0) strncpy(m_szname, sz, l);
	m_szname[l] = 0;
}

//-----------------------------------------------------------------------------
bool FEModelComponent::Init()
{
	return true;
}

//-----------------------------------------------------------------------------
//! This function checks if the 
bool FEModelComponent::IsActive()
{ 
	return m_bactive; 
}

//-----------------------------------------------------------------------------
void FEModelComponent::Activate()
{ 
	m_bactive = true; 
}

//-----------------------------------------------------------------------------
void FEModelComponent::Deactivate()
{ 
	m_bactive = false; 
}
