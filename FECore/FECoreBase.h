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



#pragma once
#include "FEParameterList.h"
#include "fecore_enum.h"
#include "FEProperty.h"
#include "ClassDescriptor.h"
#include <string>

//-----------------------------------------------------------------------------
class FECoreFactory;
class FEModel;

//-----------------------------------------------------------------------------
//! Base class for most classes in FECore library and the base class for all 
//! classes that can be registered with the framework.
class FECORE_API FECoreBase : public FEParamContainer
{
public:
	//! constructor
	FECoreBase(FEModel* fem);

	//! destructor
	virtual ~FECoreBase();

	//! set class name
	void SetName(const std::string& name);

	//! get the class' name
	const std::string& GetName() const;

	//! data serialization
	void Serialize(DumpStream& ar) override;

	//! Initialization
	virtual bool Init();

	//! validates all properties and parameters
	bool Validate() override;

	//! call this after the parameters are changed
	virtual bool UpdateParams();

public:
	//! return the super class id
	SUPER_CLASS_ID GetSuperClassID();

	//! return a (unique) string describing the type of this class
	//! This string is used in object creation
	const char* GetTypeStr();

	// Build the class' parameter and property list
	bool BuildClass();

public:
	//! number of parameters
	int Parameters() const;

	//! return a parameter
	virtual FEParam* FindParameter(const ParamString& s) override;

	//! return the property (or this) that owns a parameter
	FECoreBase* FindParameterOwner(void* pd);

public: // interface for getting/setting properties

	//! get the number of properties
	int Properties();

	//! Set a property
	bool SetProperty(int nid, FECoreBase* pm);

	//! Set a property via name
	bool SetProperty(const char* sz, FECoreBase* pm);

	//! return a property
	virtual FECoreBase* GetProperty(int i);

	//! find a property index ( returns <0 for error)
	int FindPropertyIndex(const char* szname);

	//! return a property (class)
	FEProperty* FindProperty(const char* sz, bool searchChildren = false);

	//! return a property from a paramstring
	FECoreBase* GetProperty(const ParamString& prop);

	//! return the number of properties defined
	int PropertyClasses() const;

	//! return a property
	FEProperty* PropertyClass(int i);

public:
	//! Get the parent of this object (zero if none)
	FECoreBase* GetParent();

	//! Get the ancestor of this class (this if none)
	FECoreBase* GetAncestor();

	//! Set the parent of this class
	void SetParent(FECoreBase* parent);

	//! return the component ID
	int GetID() const;

	//! set the component ID
	void SetID(int nid);

	//! Get the FE model
	FEModel* GetFEModel() const;

	//! set the FEModel of this class (use with caution!)
	void SetFEModel(FEModel* fem);

	static void SaveClass(DumpStream& ar, FECoreBase* p);
	static FECoreBase* LoadClass(DumpStream& ar, FECoreBase* p);

	// set parameters through a class descriptor
	bool SetParameters(const FEClassDescriptor& cd);
	bool SetParameters(const FEClassDescriptor::ClassVariable& cv);

public:
	//! Add a property
	//! Call this in the constructor of derived classes to 
	//! build the property list
	void AddProperty(FEProperty* pp, const char* sz, unsigned int flags = FEProperty::Required);

	void RemoveProperty(int i);

	void ClearProperties();

public:
	template <class T> T* ExtractProperty(bool extractSelf = true);

protected:
	// Set the factory class
	void SetFactoryClass(const FECoreFactory* fac);

public:
	const FECoreFactory* GetFactoryClass() const;

private:
	std::string		m_name;			//!< user defined name of component
	FECoreBase*		m_pParent;		//!< pointer to "parent" object (if any) (NOTE: only used by materials)
	FEModel*		m_fem;			//!< the model this class belongs to

	vector<FEProperty*>	m_Prop;		//!< list of properties

private:
	int		m_nID;			//!< component ID

	const FECoreFactory*	m_fac;	//!< factory class that instantiated this class

	friend class FECoreFactory;
};

// include template property definitions
#include "FEPropertyT.h"

template <class T>	FEProperty* AddClassProperty(FECoreBase* pc, T* pp, const char* sz)
{
	FEFixedPropertyT<T>* prop = new FEFixedPropertyT<T>(pp);
	prop->SetDefaultType(sz);
	pc->AddProperty(prop, sz, FEProperty::Fixed);
	return prop;
}

template <class T> FEProperty* AddClassProperty(FECoreBase* pc, T** pp, const char* sz, unsigned int flags = FEProperty::Required)
{
	FEPropertyT<T>* prop = new FEPropertyT<T>(pp);
	if (prop->GetSuperClassID() == FECLASS_ID) prop->SetDefaultType(sz);
	pc->AddProperty(prop, sz, flags);
	return prop;
}

template <class T>	FEProperty* AddClassProperty(FECoreBase* pc, std::vector<T*>* pp, const char* sz, unsigned int flags = FEProperty::Required)
{
	FEVecPropertyT<T>* prop = new FEVecPropertyT<T>(pp);
	if (prop->GetSuperClassID() == FECLASS_ID) prop->SetDefaultType(sz);
	pc->AddProperty(prop, sz, flags);
	return prop;
}

#define ADD_PROPERTY(theProp, ...) AddClassProperty(this, &theProp, __VA_ARGS__)

#define FECORE_SUPER_CLASS(a) public: static SUPER_CLASS_ID superClassID() { return a; }
//#define REGISTER_SUPER_CLASS(theClass, a) SUPER_CLASS_ID theClass::superClassID() { return a;}

template <class T> T* FECoreBase::ExtractProperty(bool extractSelf)
{
	if (extractSelf)
	{
		if (dynamic_cast<T*>(this)) return dynamic_cast<T*>(this);
	}

	int NC = Properties();
	for (int i = 0; i < NC; i++)
	{
		FECoreBase* pci = GetProperty(i);
		if (pci)
		{
			T* pc = pci->ExtractProperty<T>();
			if (pc) return pc;
		}
	}
	return nullptr;
}
