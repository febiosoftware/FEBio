#pragma once
#include "fecore_enum.h"
#include "FECoreBase.h"

//-----------------------------------------------------------------------------
//! Forward declaration of the FEModel class. All classes that register
//! with the framework take a pointer to FEModel as their constructor parameter.
class FEModel;

//-----------------------------------------------------------------------------
//! The factory class contains the mechanism for instantiating a class.
class FECORE_API FECoreFactory : public FEParamContainer
{
public:
	//! constructor
	FECoreFactory(SUPER_CLASS_ID scid, const char* sztype, int nspec = -1);

	//! virtual constructor
	virtual ~FECoreFactory();

	//! This is the function that the kernel will use to intantiate an object
	void* CreateInstance(FEModel* pfem);

public:
	// return the type string identifier
	const char* GetTypeStr() const { return m_sztype; }

	//! return the super-class ID
	SUPER_CLASS_ID GetSuperClassID() const { return m_scid; }

	//! return the module name
	unsigned int GetModuleID() const { return m_module; }

	//! set the module name
	void SetModuleID(unsigned int nid);

	//! Get the spec number
	int GetSpecID() const { return m_spec; }
	
public:
	//! derived classes implement this to create an instance of a class
	virtual void* Create(FEModel*) = 0;

private:
	const char*		m_sztype;	//!< class type string
	int				m_spec;		//!< The max spec number for which this feature is defined (-1 is don't care)
	unsigned int	m_module;	//!< ID of module this class belongs to
	SUPER_CLASS_ID	m_scid;		//!< the super-class ID
};

//-----------------------------------------------------------------------------
//! Forward declarations of classes used by the domain factory
class FEDomain;
class FEMesh;
class FEMaterial;
class FEModel;

//-----------------------------------------------------------------------------
//! Creation of domains are a little more elaborate and deviate from the usual
//! factory methods.
class FECORE_API FEDomainFactory
{
public:
	FEDomainFactory(){}
	virtual ~FEDomainFactory(){}

	virtual FEDomain* CreateDomain(const FE_Element_Spec& spec, FEMesh* pm, FEMaterial* pmat) = 0;
};
