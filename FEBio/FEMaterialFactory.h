#pragma once

#include <list>
using namespace std;

// forward declaration of the material class
class FEMaterial;

// forward declaration of the material registration class
class FERegisterMaterial;

//-----------------------------------------------------------------------------
//! Factory class for dynamic creation of materials

//! The FEMaterialFactory class allows the framework to create a class, based
//! on a text string. In order for this to work, each material derived from
//! FEMaterial should have the DECLARE_REGISTERED macro in its declaration and
//! should also have the REGISTER_MATERIAL declared somewhere. These macros
//! hide the details on how the material registration procedure works. 
//! Note that this class is a singleton, i.e. only one instance can be created.
//! Since the interface is static, you should never create an object of this
//! class explicitly. If you do need a pointer to the actual instance, call
//! the GetInstance() member.

class FEMaterialFactory
{
protected:
	typedef FEMaterial* (*MAT_CREATE_FUNC)();	// declares a material create function type

	enum {MAX_MATERIAL_NAME = 64};	// max size of material name

	struct REGISTERED_MATERIAL
	{
		char	m_szname[MAX_MATERIAL_NAME];	// name associated with material
		MAT_CREATE_FUNC	m_fnc;					// function to call to create this material
	};

	class Exception {};

public:
	//! Get the one and only instance of FEMaterialFactory
	static FEMaterialFactory* GetInstance();

	//! Register a material
	static void RegisterMaterial(const char* sz, MAT_CREATE_FUNC fnc);

	//! Create a material class based on it's registered name
	static FEMaterial* CreateMaterial(const char* szmat);

protected:
	list<REGISTERED_MATERIAL>	m_ml;	//!< the list of registered materials

	//! The one and only instance
	static FEMaterialFactory*	m_pMF;

private:
	//! default constructor
	/*! The constructor is declared as private so
	    that the user cannot instantiate it       */
	FEMaterialFactory(void) {}

	friend class FERegisterMaterial;
};

//-----------------------------------------------------------------------------
//! Registers materials with the FEMaterialFactory

//! This class helps with the registration of a material with the FEMaterialFactory
class FERegisterMaterial
{
public:
	FERegisterMaterial(const char* sz, FEMaterialFactory::MAT_CREATE_FUNC fnc)
	{
		FEMaterialFactory::RegisterMaterial(sz, fnc);
	}
};

// The following macros should be used to register a material with the material factory.
// To register a material, take the following steps.
// 1) add the DECLARE_REGISTERED macro to the material class declaration
// 2) add the REGISTER_MATERIAL macro to the material class definition

// The DECLARE_REGISTERED macro sets up the mechanism to do the material registration
#define DECLARE_REGISTERED(theClass) \
public: \
	const char* GetTypeString(); \
	static FEMaterial* CreateMaterial() {return new theClass; }\
	static FERegisterMaterial	m_##theClass##_rm;

// the REGISTER_MATERIAL does the actual material registration
#define REGISTER_MATERIAL(theClass, theName) \
	FERegisterMaterial theClass::m_##theClass##_rm(theName, theClass::CreateMaterial); \
	const char* theClass::GetTypeString() { return theName; }

