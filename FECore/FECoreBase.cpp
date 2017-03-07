#include "stdafx.h"
#include "FECoreBase.h"
#include "DumpStream.h"
#include "FECoreKernel.h"

//-----------------------------------------------------------------------------
//! The constructor takes one argument, namely the SUPER_CLASS_ID which
//! defines the type of class this is. (The SUPER_CLASS_ID was introduced to
//! eliminate a lot of akward dynamic_casts.)
FECoreBase::FECoreBase(SUPER_CLASS_ID sid) : m_sid(sid) 
{ 
	m_sztype = 0; 
	m_szname[0] = 0;
	m_pParent = 0;
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
	// do base class first
	FEParamContainer::Serialize(ar);

	// serialize name
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

	// serialize all the properties
	int NP = m_Prop.size();
	for (int i = 0; i<NP; ++i)
	{
		FEProperty* pmat = m_Prop[i];
		pmat->SetParent(this);
		pmat->Serialize(ar);
	}
}

//-----------------------------------------------------------------------------
bool FECoreBase::Validate()
{
	// call base class first
	if (FEParamContainer::Validate() == false) return false;

	// check properties
	const int nprop = (int)m_Prop.size();
	for (int i = 0; i<nprop; ++i)
	{
		FEProperty* pi = m_Prop[i];
		if (pi)
		{
			if (pi->Validate() == false) return false;
		}
	}

	return true;
}

//-----------------------------------------------------------------------------
bool FECoreBase::Init()
{
	// check the parameter ranges
	if (Validate() == false) return false;

	// initialize properties
	const int nprop = (int)m_Prop.size();
	for (int i = 0; i<nprop; ++i)
	{
		FEProperty* pi = m_Prop[i];
		if (pi)
		{
			if (pi->Init() == false) return false;
		}
		else return fecore_error("A nullptr was set for property i");
	}
	return true;
}

//-----------------------------------------------------------------------------
void FECoreBase::AddProperty(FEProperty* pp, const char* sz, unsigned int flags)
{
	pp->SetName(sz);
	pp->m_brequired = ((flags & FEProperty::Required) != 0);
	pp->m_bvalue    = ((flags & FEProperty::ValueProperty) != 0);
	m_Prop.push_back(pp);
}

//-----------------------------------------------------------------------------
int FECoreBase::Properties()
{
	int N = (int)m_Prop.size();
	int n = 0;
	for (int i = 0; i<N; ++i) n += m_Prop[i]->size();
	return n;
}

//-----------------------------------------------------------------------------
int FECoreBase::FindPropertyIndex(const char* sz)
{
	int NP = (int)m_Prop.size();
	for (int i = 0; i<NP; ++i)
	{
		const FEProperty* pm = m_Prop[i];
		if (pm && (strcmp(pm->GetName(), sz) == 0)) return i;
	}
	return -1;
}

//-----------------------------------------------------------------------------
FEProperty* FECoreBase::FindProperty(const char* sz)
{
	int NP = (int)m_Prop.size();
	for (int i = 0; i<NP; ++i)
	{
		FEProperty* pm = m_Prop[i];
		if (pm && (strcmp(pm->GetName(), sz) == 0)) return pm;
	}
	return nullptr;
}

//-----------------------------------------------------------------------------
FECoreBase* FECoreBase::GetProperty(int n)
{
	int N = (int)m_Prop.size();
	int m = 0;
	for (int i = 0; i<N; ++i)
	{
		FEProperty* pm = m_Prop[i];
		int l = pm->size();
		if (m + l > n) return pm->get(n - m);
		m += l;
	}
	return 0;
}

//-----------------------------------------------------------------------------
bool FECoreBase::SetProperty(int i, FECoreBase* pb)
{
	FEProperty* pm = m_Prop[i];
	if (pm->IsType(pb))
	{
		pm->SetProperty(pb);
		if (pb) pb->SetParent(this);
		return true;
	}
	return false;
}

//-----------------------------------------------------------------------------
FEParam* FECoreBase::GetParameter(const ParamString& s)
{
	// first search the parameter list
	FEParam* p = FEParamContainer::GetParameter(s);
	if (p) return p;

	// next, let's try the property list
	int NP = (int)m_Prop.size();
	for (int i = 0; i<NP; ++i)
	{
		// get the property
		FEProperty* mp = m_Prop[i];

		// see if matches
		if (s == mp->GetName())
		{
			if (mp->IsArray())
			{
				// get the number of items in this property
				int nsize = mp->size();
				int nid = s.Index();
				if ((nid >= 0) && (nid < nsize))
				{
					return mp->get(nid)->GetParameter(s.next());
				}
				else if (s.IndexString())
				{
					FECoreBase* c = mp->get(s.IndexString());
					if (c) return c->GetParameter(s.next());
				}
			}
			else
			{
				return mp->get(0)->GetParameter(s.next());
			}
		}
	}

	return 0;
}
