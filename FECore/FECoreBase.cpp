#include "stdafx.h"
#include "FECoreBase.h"
#include "DumpStream.h"
#include "FECoreKernel.h"
#include "log.h"

//-----------------------------------------------------------------------------
//! The constructor takes one argument, namely the SUPER_CLASS_ID which
//! defines the type of class this is. (The SUPER_CLASS_ID was introduced to
//! eliminate a lot of akward dynamic_casts.)
FECoreBase::FECoreBase(FEModel* fem) : m_fem(fem)
{ 
	m_nID = -1;
	m_pParent = 0;
	m_fac = nullptr;
}

//-----------------------------------------------------------------------------
//! destructor does nothing for now.
FECoreBase::~FECoreBase(){}

//-----------------------------------------------------------------------------
//! return the super class id
SUPER_CLASS_ID FECoreBase::GetSuperClassID() { return (m_fac ? m_fac->GetSuperClassID() : FEINVALID_ID); }

//-----------------------------------------------------------------------------
//! return a (unique) string describing the type of this class
//! This string is used in object creation
const char* FECoreBase::GetTypeStr() { return (m_fac ? m_fac->GetTypeStr() : nullptr); }

//-----------------------------------------------------------------------------
//! Set the factory class
void FECoreBase::SetFactoryClass(FECoreFactory* fac)
{
	m_fac = fac;
}

//-----------------------------------------------------------------------------
FECoreFactory* FECoreBase::GetFactoryClass()
{
	return m_fac;
}

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
//! Get the parent of this object (zero if none)
FECoreBase* FECoreBase::GetParent()
{ 
	return m_pParent; 
}

//-----------------------------------------------------------------------------
FECoreBase* FECoreBase::GetAncestor()
{
	FECoreBase* mp = GetParent(); 
	return (mp ? mp->GetAncestor() : this); 
}

//-----------------------------------------------------------------------------
//! Set the parent of this class
void FECoreBase::SetParent(FECoreBase* parent) { m_pParent = parent; }

//-----------------------------------------------------------------------------
//! return the component ID
int FECoreBase::GetID() const { return m_nID; }

//-----------------------------------------------------------------------------
//! set the component ID
void FECoreBase::SetID(int nid) { m_nID = nid; }

//-----------------------------------------------------------------------------
//! Get the FE model
FEModel* FECoreBase::GetFEModel() const { return m_fem; }

//-----------------------------------------------------------------------------
void FECoreBase::Serialize(DumpStream& ar)
{
	// do base class first
	FEParamContainer::Serialize(ar);

	// serialize name
	if (ar.IsShallow() == false)
	{
		ar & m_name;
		ar & m_nID;
	}

	// serialize all the properties
	int NP = (int)m_Prop.size();
	for (int i = 0; i<NP; ++i)
	{
		FEProperty* prop = m_Prop[i];
		prop->SetParent(this);
		prop->Serialize(ar);
	}
}

//-----------------------------------------------------------------------------
void FECoreBase::SaveClass(DumpStream& ar, FECoreBase* a)
{
	assert(ar.IsSaving());
	int classID = 0;
	if (a == nullptr) { ar << classID; return; }
	classID = a->GetSuperClassID();
	const char* sztype = a->GetTypeStr();
	ar << classID;
	ar << sztype;
}

FECoreBase* FECoreBase::LoadClass(DumpStream& ar, FECoreBase* a)
{
	assert(ar.IsLoading());

	int classID = 0;
	ar >> classID;
	if (classID == FEINVALID_ID) return nullptr;

	char sztype[256] = { 0 };
	ar >> sztype;

	// instantiate the class
	a = fecore_new<FECoreBase>(classID, sztype, &ar.GetFEModel());
	assert(a);

	if (a == nullptr) throw DumpStream::ReadError();

	return a;
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
			if (pi->Init() == false)
			{
				feLogError("The required property \"%s\" was not defined", pi->GetName());
				return false;
			}
		}
		else
		{
			feLogError("A nullptr was set for property i");
			return false;
		}
	}
	return true;
}

//-----------------------------------------------------------------------------
void FECoreBase::AddProperty(FEProperty* pp, const char* sz, unsigned int flags)
{
	pp->SetName(sz);
	pp->SetFlags(flags);
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
//! number of parameters
int FECoreBase::Parameters() const
{
	return GetParameterList().Parameters();
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
//! return the property (or this) that owns a parameter
FECoreBase* FECoreBase::FindParameterOwner(void* pd)
{
	// see if this class is the owner of the data pointer
	FEParam* p = FindParameterFromData(pd);
	if (p) return this;

	// it's not se let's check the properties
	int NP = PropertyClasses();
	for (int i = 0; i < NP; ++i)
	{
		FEProperty* pi = PropertyClass(i);
		int n = pi->size();
		for (int j = 0; j < n; ++j)
		{
			FECoreBase* pcj = pi->get(j);
			if (pcj)
			{
				FECoreBase* pc = pcj->FindParameterOwner(pd);
				if (pc) return pc;
			}
		}
	}

	// sorry, no luck
	return nullptr;
}

//-----------------------------------------------------------------------------
//! return the number of properties defined
int FECoreBase::PropertyClasses() const 
{ 
	return (int)m_Prop.size(); 
}

//! return a property
FEProperty* FECoreBase::PropertyClass(int i)
{ 
	return m_Prop[i]; 
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
