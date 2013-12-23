#pragma once
#include "FEParameterList.h"
#include "FE_enum.h"

//-----------------------------------------------------------------------------
class FEBioFactory;

//-----------------------------------------------------------------------------
//! Base class for most classes in FECore library and the base class for all 
//! classes that can be registered with the framework.
class FECoreBase : public FEParamContainer
{
public:
	FECoreBase(SUPER_CLASS_ID sid);
	virtual ~FECoreBase();

public:
	//! return the super class id
	SUPER_CLASS_ID GetSuperClassID();

	//! return a (unique) string describing the type of this class
	//! This string is used in object creation
	const char* GetTypeStr();

public: // interface for managing attributes

	//! Set the attribute
	virtual bool SetAttribute(const char* szname, const char* szval) { return true; }

private:
	//! Set the type string (This is used by the factory methods to make sure 
	//! the class has the same type string as corresponding factory class
	void SetTypeStr(const char* sz);

private:
	SUPER_CLASS_ID	m_sid;		//!< The super-class ID
	const char*		m_sztype;	//!< the type string

	friend class FEBioFactory;
};
