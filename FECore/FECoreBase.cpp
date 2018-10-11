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
	m_nID = -1;
	m_sztype = 0;
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
void FECoreBase::SetName(const std::string& name)
{ 
	m_name = name;
}

//-----------------------------------------------------------------------------
//! Return the name
const std::string& FECoreBase::GetName() const
{ 
	return m_name; 
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
			ar << m_name;
			ar << m_nID;
		}
		else
		{
			ar >> m_name;
			ar >> m_nID;
		}
	}

	// serialize all the properties
	int NP = (int)m_Prop.size();
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
			if (pi->Init() == false) return fecore_error("The required property \"%s\" was not defined", pi->GetName());
		}
		else return fecore_error("A nullptr was set for property i");
	}
	return true;
}

//-----------------------------------------------------------------------------
void FECoreBase::AddProperty(FEProperty* pp, const char* sz, unsigned int flags)
{
	pp->SetName(sz);
	pp->m_brequired = ((flags & FEProperty::Optional) == 0);
	pp->SetParent(this);
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
FEParam* FECoreBase::FindParameter(const ParamString& s)
{
	// first search the parameter list
	FEParam* p = FEParamContainer::FindParameter(s);
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
				int index = s.Index();
				if ((index >= 0) && (index < nsize))
				{
					return mp->get(index)->FindParameter(s.next());
				}
				else
				{
					int nid = s.ID();
					if (nid != -1)
					{
						FECoreBase* pc = mp->getFromID(nid);
						if (pc) return pc->FindParameter(s.next());
					}
					else if (s.IDString())
					{
						FECoreBase* c = mp->get(s.IDString());
						if (c) return c->FindParameter(s.next());
					}
				}
			}
			else
			{
				FECoreBase* pc = mp->get(0);
				return (pc ? pc->FindParameter(s.next()) : nullptr);
			}
		}
	}

	return nullptr;
}

//-----------------------------------------------------------------------------
FECoreBase* FECoreBase::GetProperty(const ParamString& prop)
{
	int NP = (int) m_Prop.size();
	for (int i=0; i<NP; ++i)
	{	
		FEProperty* mp = m_Prop[i];

		if (prop == mp->GetName())
		{
			if (mp->IsArray())
			{
				// get the number of items in this property
				int nsize = mp->size();
				int index = prop.Index();
				if ((index >= 0) && (index < nsize))
				{
					FECoreBase* pc = mp->get(index);
					if (pc)
					{
						ParamString next = prop.next();
						if (next.count() == 0) return pc;
						else return pc->GetProperty(next);
					}
				}
				else
				{
					int nid = prop.ID();
					if (nid != -1)
					{
						FECoreBase* pc = mp->getFromID(nid);
					}
					else if (prop.IDString())
					{
						FECoreBase* pc = mp->get(prop.IDString());
						if (pc)
						{
							ParamString next = prop.next();
							if (next.count() == 0) return pc;
							else return pc->GetProperty(next);
						}
					}
				}
			}
			else
			{
				FECoreBase* pc = mp->get(0);
				ParamString next = prop.next();
				if (next.count() == 0) return pc;
				else return pc->GetProperty(next);
			}
		}
	}

	return 0;
}

//-----------------------------------------------------------------------------
bool FECoreBase::BuildClass()
{
	GetParameterList();

	for (int i = 0; i < Properties(); ++i)
	{
		FEProperty* pp = PropertyClass(i);
		int m = pp->size();
		for (int j = 0; j < m; ++j)
		{
			FECoreBase* pj = pp->get(j);
			if (pj) pj->BuildClass();
		}
	}
	return true;
}
