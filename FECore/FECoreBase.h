#pragma once
#include "FEParameterList.h"
#include "FE_enum.h"

//-----------------------------------------------------------------------------
class FECoreFactory;

//-----------------------------------------------------------------------------
//! Base class for most classes in FECore library and the base class for all 
//! classes that can be registered with the framework.
class FECoreBase : public FEParamContainer
{
public:
	//! constructor
	FECoreBase(SUPER_CLASS_ID sid);

	//! destructor
	virtual ~FECoreBase();

	//! set material name
	void SetName(const char* sz);

	//! get the material's name
	const char* GetName();

	//! data serialization
	void Serialize(DumpStream& ar);

public:
	//! return the super class id
	SUPER_CLASS_ID GetSuperClassID();

	//! return a (unique) string describing the type of this class
	//! This string is used in object creation
	const char* GetTypeStr();

public: // interface for managing attributes

	//! Set the attribute
	virtual bool SetAttribute(const char* szname, const char* szval) { return true; }

public: // interface for getting/setting properties

	//! get the number of properties
	virtual int Properties () { return 0; }

	//! get a specific property
	virtual FECoreBase* GetProperty(int i) { return 0; }

	//! find a property index ( returns <0 for error)
	virtual int FindPropertyIndex(const char* szname) { return -1; }

	//! set a property (returns false on error)
	virtual bool SetProperty(int i, FECoreBase* pm) { return false; }

private:
	//! Set the type string (This is used by the factory methods to make sure 
	//! the class has the same type string as corresponding factory class
	void SetTypeStr(const char* sz);

private:
	char			m_szname[128];	//!< user defined name of component

	SUPER_CLASS_ID	m_sid;		//!< The super-class ID
	const char*		m_sztype;	//!< the type string

	friend class FECoreFactory;
};
