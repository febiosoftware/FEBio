#pragma once
#include "FE_enum.h"
#include "FECoreBase.h"

//-----------------------------------------------------------------------------
//! Forward declaration of the FEModel class. All classes that register
//! with the framework take a pointer to FEModel as their constructor parameter.
class FEModel;

//-----------------------------------------------------------------------------
//! The factory class contains the mechanism for instantiating a class.
class FECoreFactory
{
public:
	//! constructor
	FECoreFactory(SUPER_CLASS_ID scid, const char* sztype);

	//! virtual constructor
	virtual ~FECoreFactory();

	//! This is the function that the kernel will use to intantiate an object
	void* CreateInstance(FEModel* pfem);

public:
	// return the type string identifier
	const char* GetTypeStr() const { return m_sztype; }

	//! return the super-class ID
	SUPER_CLASS_ID GetSuperClassID() const { return m_scid; }

public:
	//! derived classes implement this to create an instance of a class
	virtual void* Create(FEModel*) = 0;

private:
	const char*		m_sztype;	//!< class type string
	SUPER_CLASS_ID	m_scid;		//!< the super-class ID
};

//-----------------------------------------------------------------------------
//! Forward declarations of classes used by the domain factory
class FEDomain;
class FEMesh;
class FEMaterial;

//-----------------------------------------------------------------------------
//! Creation of domains are a little more elaborate and deviate from the usual
//! factory methods.
class FEDomainFactory
{
public:
	FEDomainFactory(){}
	virtual ~FEDomainFactory(){}

	virtual FEDomain* CreateDomain(const FE_Element_Spec& spec, FEMesh* pm, FEMaterial* pmat) = 0;
};

//-----------------------------------------------------------------------------
// factory class for linear solvers.
class LinearSolver;

class FELinearSolverFactory
{
public:
	FELinearSolverFactory(int nid) : m_nsolver_id(nid) {}
	virtual ~FELinearSolverFactory(){}

	virtual LinearSolver* Create() = 0;

	int GetID() { return m_nsolver_id; }

private:
	int	m_nsolver_id;
};
