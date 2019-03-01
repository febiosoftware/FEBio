#pragma once
#include "FEParameterList.h"
#include "fecore_enum.h"
#include "FEProperty.h"
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
	void Serialize(DumpStream& ar);

	//! Initialization
	virtual bool Init();

	//! validates all properties and parameters
	bool Validate();

public:
	//! return the super class id
	SUPER_CLASS_ID GetSuperClassID();

	//! return a (unique) string describing the type of this class
	//! This string is used in object creation
	const char* GetTypeStr();

	// Build the class' parameter and property list
	bool BuildClass();

public:
	//! return a parameter
	virtual FEParam* FindParameter(const ParamString& s);

	//! return the property (or this) that owns a parameter
	FECoreBase* FindParameterOwner(void* pd);

public: // interface for getting/setting properties

	//! get the number of properties
	int Properties();

	//! Set a property
	bool SetProperty(int nid, FECoreBase* pm);

	//! return a property
	virtual FECoreBase* GetProperty(int i);

	//! find a property index ( returns <0 for error)
	int FindPropertyIndex(const char* szname);

	//! return a property (class)
	FEProperty* FindProperty(const char* sz);

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

	static void SaveClass(DumpStream& ar, FECoreBase* p);
	static FECoreBase* LoadClass(DumpStream& ar, FECoreBase* p);

public:
	//! Add a property
	//! Call this in the constructor of derived classes to 
	//! build the property list
	void AddProperty(FEProperty* pp, const char* sz, unsigned int flags = FEProperty::Required);

public:
	template <class T> T* ExtractProperty(bool extractSelf = true);

protected:
	// Set the factory class
	void SetFactoryClass(FECoreFactory* fac);

public:
	FECoreFactory* GetFactoryClass();

private:
	std::string		m_name;			//!< user defined name of component
	FECoreBase*		m_pParent;		//!< pointer to "parent" object (if any) (NOTE: only used by materials)
	FEModel*		m_fem;			//!< the model this class belongs to

	vector<FEProperty*>	m_Prop;		//!< list of properties

private:
	int		m_nID;			//!< component ID

	FECoreFactory*	m_fac;	//!< factory class that instantiated this class

	friend class FECoreFactory;
};

// include template property definitions
#include "FEPropertyT.h"

template <class T>	void AddClassProperty(FECoreBase* pc, T** pp, const char* sz, unsigned int flags = FEProperty::Required)
{
	FEPropertyT<T>* prop = new FEPropertyT<T>(pp);
	pc->AddProperty(prop, sz, flags);
}

template <class T>	void AddClassProperty(FECoreBase* pc, std::vector<T*>* pp, const char* sz, unsigned int flags = FEProperty::Required)
{
	FEVecPropertyT<T>* prop = new FEVecPropertyT<T>(pp);
	pc->AddProperty(prop, sz, flags);
}

#define ADD_PROPERTY(theProp, ...) AddClassProperty(this, &theProp, __VA_ARGS__);

#define FECORE_SUPER_CLASS public: static SUPER_CLASS_ID classID();
#define REGISTER_SUPER_CLASS(theClass, a) SUPER_CLASS_ID theClass::classID() { return a;}

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
