#include "stdafx.h"
#include "FECoreBase.h"
#include "DumpStream.h"

//-----------------------------------------------------------------------------
//! The constructor takes one argument, namely the SUPER_CLASS_ID which
//! defines the type of class this is. (The SUPER_CLASS_ID was introduced to
//! eliminate a lot of akward dynamic_casts.)
FECoreBase::FECoreBase(SUPER_CLASS_ID sid) : m_sid(sid) 
{ 
	m_sztype = 0; 
	m_szname[0] = 0;
}

//-----------------------------------------------------------------------------
//! destructor does nothing for now.
FECoreBase::~FECoreBase(){}

//-----------------------------------------------------------------------------
//! return the super class id
SUPER_CLASS_ID FECoreBase::GetSuperClassID() { return m_sid; }

//-----------------------------------------------------------------------------
//! return a (unique) string describing the type of this class
//! This string is used in object creation
const char* FECoreBase::GetTypeStr() { return m_sztype; }

//-----------------------------------------------------------------------------
//! Set the type string (This is used by the factory methods to make sure 
//! the class has the same type string as corresponding factory class
void FECoreBase::SetTypeStr(const char* sz) { m_sztype = sz; }

//-----------------------------------------------------------------------------
//! Sets the user defined name of the component
void FECoreBase::SetName(const char* sz)
{ 
	if (sz==0) return;
	strcpy(m_szname, sz); 
}

//-----------------------------------------------------------------------------
//! Return the name
const char* FECoreBase::GetName()
{ 
	return m_szname; 
}

//-----------------------------------------------------------------------------
void FECoreBase::Serialize(DumpStream& ar)
{
	FEParamContainer::Serialize(ar);
	if (ar.IsShallow() == false)
	{
		if (ar.IsSaving())
		{
			ar << m_szname;
		}
		else
		{
			ar >> m_szname;
		}
	}
}
