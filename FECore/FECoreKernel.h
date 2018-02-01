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
class FECOREDLL_EXPORT FECoreKernel
{
public:
	// Do not call this function from a plugin as it will not return the correct
	// instance. Instead, use the FECoreKernel object that is passed in the PluginInitialize method
	static FECoreKernel& GetInstance();

	// set the instance of the kernel
	static void SetInstance(FECoreKernel* pkernel);

	// Get the logfile
	static Logfile& GetLogfile();

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

	//! find a factory class
	FECoreFactory* FindFactoryClass(int classID, const char* sztype);

	//! set the active module
	void SetActiveModule(const char* szmodule);

	//! Get the active module
	const char* GetActiveModule() const;

public:
	//! Register a new domain class
	void RegisterDomain(FEDomainFactory* pf);

	//! Create a domain of a certain type
	FEDomain* CreateDomain(const FE_Element_Spec& spec, FEMesh* pm, FEMaterial* pmat);

public: // error reporting mechanism

	// Set the error string
	void SetErrorString(const char* sz);

	// returns error string (can be null if no error reported)
	const char* GetErrorString();

public:
	//! Register a linear solver factory
	void RegisterLinearSolver(FELinearSolverFactory* pf);

	//! create a linear solver factory
	LinearSolver* CreateLinearSolver(int nsolver);

	//! Find linear solver factory
	FELinearSolverFactory* FindLinearSolverFactory(int nsolver);

public:
	//! set the default linear solver
	static void SetDefaultSolver(int nsolver) { m_ndefault_solver = nsolver; }
	static int m_ndefault_solver;

private:
	std::vector<FECoreFactory*>			m_Fac;	// list of registered factory classes
	std::vector<FEDomainFactory*>		m_Dom;	// list of domain factory classes
	std::vector<FELinearSolverFactory*> m_LS;	// list of linear solver factories

	// active module name
	const char* m_szmod;

	Logfile*	m_plog;	// keep a pointer to the logfile (used by plugins)

	char*	m_szerr;	//!< error string

private: // make singleton
	FECoreKernel();
	FECoreKernel(const FECoreKernel&){}
	void operator = (const FECoreKernel&){}

private:
	static FECoreKernel* m_pKernel;	// the one-and-only kernel object
};

//-----------------------------------------------------------------------------
// helper function for reporting errors
bool fecore_error(const char* szerr, ...);
const char* fecore_get_error_string();

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
// Register a deprecated class using default creation parameters
#define REGISTER_FECORE_CLASS_OBSOLETE(theClass, theSID, theName) \
	static FERegisterClass_T<theClass> _##theClass##_old_rc(theSID, theName);

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

//=============================================================================
// TODO: Move all this stuff to sdk.h

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
