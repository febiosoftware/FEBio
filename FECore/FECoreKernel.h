#pragma once
#include "FECoreFactory.h"
#include <vector>
#include <string.h>
#include <stdio.h>
#include "version.h"

//-----------------------------------------------------------------------------
// Forward declaration of model class
class FEModel;

//-----------------------------------------------------------------------------
// forward declaration of the log file
class Logfile;

//-----------------------------------------------------------------------------
//! This is the FECore kernel class that manages the interactions between the 
//! different modules. In particular, it manages the factory classes
//! which are responsible for the creation of different classes that are registered
//! with the kernel.
class FECoreKernel
{
public:
	// Do not call this function from a plugin as it will not return the correct
	// instance. Instead, use the FECoreKernel object that is passed in the PluginInitialize method
	static FECoreKernel& GetInstance();

	// set the instance of the kernel
	static void SetInstance(FECoreKernel* pkernel);

	// Get the logfile
	Logfile& GetLogfile();

public:
	//! Register a class with the framework
	void RegisterClass(FECoreFactory* ptf);

	//! Create a specific class
	void* Create(SUPER_CLASS_ID, const char* sztag, FEModel* pfem);

	//! count the number of registered classes with a given super-class id
	int Count(SUPER_CLASS_ID sid);

	//! List the registered classes with a given super-class id
	void List(SUPER_CLASS_ID sid);

	//! Get the number of registered factory classes
	int FactoryClasses();

	//! return a factory class
	const FECoreFactory* GetFactoryClass(int i);

public:
	//! Register a new domain class
	void RegisterDomain(FEDomainFactory* pf);

	//! Create a domain of a certain type
	FEDomain* CreateDomain(const FE_Element_Spec& spec, FEMesh* pm, FEMaterial* pmat);

public:
	//! Register a linear solver factory
	void RegisterLinearSolver(FELinearSolverFactory* pf);

	//! create a linear solver factory
	LinearSolver* CreateLinearSolver(int nsolver);

private:
	std::vector<FECoreFactory*>			m_Fac;	// list of registered factory classes
	std::vector<FEDomainFactory*>		m_Dom;	// list of domain factory classes
	std::vector<FELinearSolverFactory*> m_LS;	// list of linear solver factories

	Logfile*	m_plog;	// keep a pointer to the logfile (used by plugins)

private: // make singleton
	FECoreKernel();
	FECoreKernel(const FECoreKernel&){}
	void operator = (const FECoreKernel&){}

private:
	static FECoreKernel* m_pKernel;	// the one-and-only kernel object
};

//-----------------------------------------------------------------------------
//! This class helps with the registration of a class with the framework
template <typename T> class FERegisterClass_T : public FECoreFactory
{
public:
	FERegisterClass_T(SUPER_CLASS_ID sid, const char* sz) : FECoreFactory(sid, sz)
	{
		FECoreKernel& fecore = FECoreKernel::GetInstance();
		fecore.RegisterClass(this);
	}
	void* Create(FEModel* pfem) { return new T(pfem); }
};

//-----------------------------------------------------------------------------
// Register a class using default creation parameters
#define REGISTER_FECORE_CLASS(theClass, theSID, theName) \
	static FERegisterClass_T<theClass> _##theClass##_rc(theSID, theName);

//-----------------------------------------------------------------------------
// version for classes that require template arguments
#define REGISTER_FECORE_CLASS_T(theClass, theSID, theArg, theName) \
	static FERegisterClass_T<theClass<theArg> > _##theClass##theArg##_rc(theSID, theName);

//-----------------------------------------------------------------------------
template <typename TBase> inline TBase* fecore_new(SUPER_CLASS_ID sid, const char* sztype, FEModel* pfem)
{
	FECoreKernel& fecore = FECoreKernel::GetInstance();
	return static_cast<TBase*>(fecore.Create(sid, sztype, pfem));
}

//-----------------------------------------------------------------------------
// Template class for factory classes for plugins
template <typename T, SUPER_CLASS_ID sid> class FEPluginFactory_T : public FECoreFactory
{
public:
	FEPluginFactory_T(const char* sz) : FECoreFactory(sid, sz){}
	void* Create(FEModel* pfem) { return new T(pfem); }
};

//------------------------------------------------------------------------------
// This is for functions exported from a plugin
#ifdef WIN32
#define FECORE_EXPORT extern "C" __declspec(dllexport)
#else
#define FECORE_EXPORT extern "C"
#endif
