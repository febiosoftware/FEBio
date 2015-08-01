#include "stdafx.h"
#include "FEModelComponent.h"
#include <string.h>

//-----------------------------------------------------------------------------
//! The constructor takes two arguments: the SUPER_CLASS_ID which defines the 
//! type of the model component and a pointer to the FEModel object this component
//! belongs to.
FEModelComponent::FEModelComponent(SUPER_CLASS_ID sid, FEModel* pfem) : FECoreBase(sid)
{
	// assign a class ID
	static int nid = 1;
	m_nClassID = nid++;

	// the ID can be used by derived class to define a identifier for derived classes
	// This value needs to be set in the constructor
	m_nID = 0;

	// initialize parameters
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
FEModel* FEModelComponent::GetFEModel() const
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
int FEModelComponent::GetID() const
{
	return m_nID;
}

//-----------------------------------------------------------------------------
void FEModelComponent::SetID(int n)
{
	m_nID = n;
}

//-----------------------------------------------------------------------------
int FEModelComponent::GetClassID() const
{
	return m_nClassID;
}

//-----------------------------------------------------------------------------
void FEModelComponent::SetClassID(int n)
{
	m_nClassID = n;
}

//-----------------------------------------------------------------------------
bool FEModelComponent::Init()
{
	return true;
}

//-----------------------------------------------------------------------------
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

//-----------------------------------------------------------------------------
void FEModelComponent::Serialize(DumpFile& ar)
{
	FECoreBase::Serialize(ar);
	if (ar.IsSaving())
	{
		ar << m_nID;
		ar << m_nClassID;
		ar << m_bactive;
		ar << m_szname;
	}
	else
	{
		char szname[256];
		ar >> m_nID;
		ar >> m_nClassID;
		ar >> m_bactive;
		ar >> szname;
		SetName(szname);
	}
}
