/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
the City of New York, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/



#include "stdafx.h"
#include "FECoreBase.h"
#include "DumpStream.h"
#include "FECoreKernel.h"
#include "FEModelParam.h"
#include "log.h"
#include <sstream>

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
FECoreBase::~FECoreBase()
{
	ClearProperties();
}

//-----------------------------------------------------------------------------
//! return the super class id
SUPER_CLASS_ID FECoreBase::GetSuperClassID() { return (m_fac ? m_fac->GetSuperClassID() : FEINVALID_ID); }

//-----------------------------------------------------------------------------
//! return a (unique) string describing the type of this class
//! This string is used in object creation
const char* FECoreBase::GetTypeStr() { return (m_fac ? m_fac->GetTypeStr() : nullptr); }

//-----------------------------------------------------------------------------
//! Set the factory class
void FECoreBase::SetFactoryClass(const FECoreFactory* fac)
{
	m_fac = fac;
}

//-----------------------------------------------------------------------------
const FECoreFactory* FECoreBase::GetFactoryClass() const
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
void FECoreBase::SetFEModel(FEModel* fem) { m_fem = fem; }

//-----------------------------------------------------------------------------
void setParamValue(FEParam& pi, const std::string& val)
{
	if (val.empty()) return;
	const char* sz = val.c_str();
	switch (pi.type())
	{
	case FE_PARAM_INT: pi.value<int>() = atoi(sz); break;
	case FE_PARAM_BOOL: pi.value<bool>() = (atoi(sz) == 0 ? false : true); break;
	case FE_PARAM_DOUBLE: pi.value<double>() = atof(sz); break;
	default:
		assert(false);
	}
}

//-----------------------------------------------------------------------------
// set parameters through a class descriptor
bool FECoreBase::SetParameters(const FEClassDescriptor& cd)
{
	const FEClassDescriptor::ClassVariable* root = cd.Root();
	return SetParameters(*cd.Root());
}

//-----------------------------------------------------------------------------
// set parameters through a class descriptor
bool FECoreBase::SetParameters(const FEClassDescriptor::ClassVariable& cv)
{
	FEParameterList& PL = GetParameterList();
	for (int i=0; i<cv.Count(); ++i)
	{
		// get the next variable
		const FEClassDescriptor::Variable* vari = cv.GetVariable(i);

		// see if this parameter is defined
		FEParam* pi = PL.FindFromName(vari->m_name.c_str());
		if (pi)
		{
			const FEClassDescriptor::SimpleVariable* vi = dynamic_cast<const FEClassDescriptor::SimpleVariable*>(vari);
			assert(vi);
			if (vi == nullptr) return false;

			// set the value
			setParamValue(*pi, vi->m_val);
		}
		else
		{
			// could be a property
			const FEClassDescriptor::ClassVariable* ci = dynamic_cast<const FEClassDescriptor::ClassVariable*>(vari);
			assert(ci);

			// find the property
			FEProperty* prop = FindProperty(ci->m_name.c_str()); assert(prop);
			if (prop == nullptr) return false;

			// allocate a new child class
			FECoreBase* pc = fecore_new<FECoreBase>(prop->GetSuperClassID(), ci->m_type.c_str(), GetFEModel()); assert(pc);
			if (pc == nullptr) return false;

			// assign the property
			prop->SetProperty(pc);

			// set the property's parameters
			pc->SetParameters(*ci);
		}
	}

	return true;
}

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

	if (ar.IsShallow() == false) ar & m_pParent;

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
	assert(classID != FEINVALID_ID);
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
	// validate parameters
	FEParameterList& pl = GetParameterList();
	int N = pl.Parameters();
	list<FEParam>::iterator pi = pl.first();
	for (int i = 0; i < N; ++i, pi++)
	{
		FEParam& p = *pi;
		if (p.is_valid() == false)
		{
			stringstream ss;
			ss << GetName() << "." << p.name();
			string paramName = ss.str();
			feLogError("Invalid value for parameter: %s", paramName.c_str());
			return false;
		}
	}

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
	// call init on model parameters
	FEParameterList& PL = GetParameterList();
	FEParamIterator it = PL.first();
	for (int i = 0; i < PL.Parameters(); ++i, ++it)
	{
		FEParam& pi = *it;
		if (pi.type() == FE_PARAM_DOUBLE_MAPPED)
		{
			for (int j = 0; j < pi.dim(); ++j)
			{
				FEParamDouble& pd = pi.value<FEParamDouble>(j);
				if (pd.Init() == false)
				{
					feLogError("Failed to initialize parameter %s", pi.name());
					return false;
				}
			}
		}
		else if (pi.type() == FE_PARAM_VEC3D_MAPPED)
		{
			for (int j = 0; j < pi.dim(); ++j)
			{
				FEParamVec3& pd = pi.value<FEParamVec3>(j);
				if (pd.Init() == false)
				{
					feLogError("Failed to initialize parameter %s", pi.name());
					return false;
				}
			}
		}
		else if (pi.type() == FE_PARAM_MAT3D_MAPPED)
		{
			for (int j = 0; j < pi.dim(); ++j)
			{
				FEParamMat3d& pd = pi.value<FEParamMat3d>(j);
				if (pd.Init() == false)
				{
					feLogError("Failed to initialize parameter %s", pi.name());
					return false;
				}
			}
		}
	}

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
				feLogError("The property \"%s\" failed to initialize or was not defined.", pi->GetName());
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
	pp->SetLongName(sz);
	pp->SetFlags(flags);
	pp->SetParent(this);
	m_Prop.push_back(pp);
}

//-----------------------------------------------------------------------------
void FECoreBase::RemoveProperty(int i)
{
	m_Prop[i] = nullptr;
}

//-----------------------------------------------------------------------------
void FECoreBase::ClearProperties()
{
	for (int i = 0; i < m_Prop.size(); ++i)
	{
		delete m_Prop[i];
	}
	m_Prop.clear();
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
FEProperty* FECoreBase::FindProperty(const char* sz, bool searchChildren)
{
	// first, search the class' properties
	int NP = (int)m_Prop.size();
	for (int i = 0; i<NP; ++i)
	{
		FEProperty* pm = m_Prop[i];
		if (pm && (strcmp(pm->GetName(), sz) == 0)) return pm;
	}

	// the property, wasn't found so look into the properties' properties
	if (searchChildren)
	{
		for (int i = 0; i < NP; ++i)
		{
			FEProperty* pm = m_Prop[i];
			if (pm)
			{
				int m = pm->size();
				for (int j = 0; j < m; ++j)
				{
					FECoreBase* pcj = pm->get(j);
					if (pcj)
					{
						// Note: we don't search children's children!
						FEProperty* pj = pcj->FindProperty(sz);
						if (pj) return pj;
					}
				}
			}
		}
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
//! Set a property via name
bool FECoreBase::SetProperty(const char* sz, FECoreBase* pb)
{
	FEProperty* prop = FindProperty(sz);
	if (prop == nullptr) return false;

	if (prop->IsType(pb))
	{
		prop->SetProperty(pb);
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

	for (int i = 0; i < PropertyClasses(); ++i)
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

//-----------------------------------------------------------------------------
bool FECoreBase::UpdateParams()
{
	return true;
}
