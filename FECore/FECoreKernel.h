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
#include "FECoreFactory.h"
#include "ClassDescriptor.h"
#include <vector>
#include <map>
#include <string.h>
#include <stdio.h>
#include "version.h"

//-----------------------------------------------------------------------------
// Forward declarations
class FEModel;
class LinearSolver;
class FEModule;


//-----------------------------------------------------------------------------
// base class for handling create events.
class FECreateHandler{
public:
	FECreateHandler() { m_moduleId = -1; }
	virtual ~FECreateHandler() {}
	virtual void handle(FECoreBase*) = 0;

	int GetModuleID() const { return m_moduleId; }
	void SetModuleID(int n) { m_moduleId = n; }

private:
	int	m_moduleId;
};

//-----------------------------------------------------------------------------
//! This is the FECore kernel class that manages the interactions between the 
//! different modules. In particular, it manages the factory classes
//! which are responsible for the creation of different classes that are registered
//! with the kernel.
class FECORE_API FECoreKernel
{
public:
	// Do not call this function from a plugin as it will not return the correct
	// instance. Instead, use the FECoreKernel object that is passed in the PluginInitialize method
	static FECoreKernel& GetInstance();

	// set the instance of the kernel
	static void SetInstance(FECoreKernel* pkernel);

public:
	static const char* SuperClassString(unsigned int sid);

	static std::map<unsigned int, const char*>	GetSuperClassMap();

public:
	//! Register a class with the framework
	void RegisterFactory(FECoreFactory* ptf);

	//! Create a specific using a superclass ID and an alias
	FECoreBase* Create(int superClassID, const char* szalias, FEModel* pfem);

	//! Create a class from its base class name and type string
	FECoreBase* Create(const char* baseClassName, const char* typeStr, FEModel* pfem);

	//! Create a specific class
	FECoreBase* CreateClass(const char* szclassName, FEModel* fem);

	//! Create a class from a class descriptor
	FECoreBase* Create(int superClassID, FEModel* pfem, const FEClassDescriptor& cd);

	//! count the number of registered classes with a given super-class id
	int Count(SUPER_CLASS_ID sid);

	//! List the registered classes with a given super-class id
	void List(SUPER_CLASS_ID sid);

	//! Get the number of registered factory classes
	int FactoryClasses();

	//! return a factory class
	const FECoreFactory* GetFactoryClass(int i);

	//! return a factory class
	const FECoreFactory* GetFactoryClass(int superClassID, int i);

	//! Get the index of a class factory (NOTE: this is a slow function!)
	int GetFactoryIndex(int superClassId, const char* sztype);

	//! find a factory class
	FECoreFactory* FindFactoryClass(int classID, const char* sztype);

	//! remove a factory class
	bool UnregisterFactory(FECoreFactory* ptf);

	//! unregister factories from allocator
	void UnregisterFactories(int alloc_id);

	//! set the current allocator ID
	void SetAllocatorID(int alloc_id);

	//! generate a allocator ID
	int GenerateAllocatorID();

	FECoreBase* CreateInstance(const FECoreFactory* fac, FEModel* fem);

	bool IsModuleActive(int moduleID);

public: // Modules

	//! count modules
	int Modules() const;

	//! Create a module (also makes it the active module)
	bool CreateModule(const char* szmodule, const char* description = nullptr);
	bool CreateModule(FEModule* pmodule, const char* szmodule, const char* description = nullptr);

	//! set the active module
	bool SetActiveModule(const char* szmodule);
	bool SetActiveModule(int moduleId);

	// return the active module's ID
	int GetActiveModuleID();
	FEModule* GetActiveModule();

	//! add module dependency to the active module
	bool AddModuleDependency(const char* szdep);

	//! Get a module's name
	const char* GetModuleName(int i) const;
	const char* GetModuleNameFromId(int id) const;
	const char* GetModuleDescription(int i) const;
	int GetModuleStatus(int i) const;

	//! Get a module's dependencies
	vector<int> GetModuleDependencies(int i) const;

	//! set the spec ID. Features with a matching spec ID will be preferred
	//! set spec ID to -1 to stop caring
	void SetSpecID(int nspec);

public:
	//! Register a new domain class
	void RegisterDomain(FEDomainFactory* pf, bool pushFront = false);

	//! Create a domain of a certain type (this uses the domain factories)
	FEDomain* CreateDomain(const FE_Element_Spec& spec, FEMesh* pm, FEMaterial* pmat);

	//! Create a domain of a certain type (this does not use the domain factories)
	FEDomain* CreateDomainExplicit(int superClass, const char* sztype, FEModel* fem);

public:
	//! set the default linear solver
	FECoreFactory* SetDefaultSolverType(const char* sztype);

	void SetDefaultSolver(FEClassDescriptor* linsolve);

	//! get the linear solver type
	const char* GetLinearSolverType() const;

	//! Get a linear solver
	LinearSolver* CreateDefaultLinearSolver(FEModel* fem);

public:
	void OnCreateEvent(FECreateHandler* pf);

	void BlockEvents(bool b);

private:
	std::vector<FECoreFactory*>			m_Fac;	// list of registered factory classes
	std::vector<FEDomainFactory*>		m_Dom;	// list of domain factory classes

	std::vector<FECreateHandler*>		m_createHandlers;
	bool								m_blockEvents;

	std::map<unsigned int, const char*>	m_sidMap;	// super class ID map

	std::string			m_default_solver_type;	// default linear solver
	FEClassDescriptor*	m_default_solver;

	// module list
	vector<FEModule*>	m_modules;
	int					m_activeModule;

	int				m_nspec;

	int		m_alloc_id;			//!< current allocator ID
	int		m_next_alloc_id;	//!< next allocator ID

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
	FERegisterClass_T(SUPER_CLASS_ID sid, const char* szclass, const char* szbase, const char* szalias, int spec = -1) : FECoreFactory(sid, szclass, szbase, szalias, spec)
	{
		FECoreKernel& fecore = FECoreKernel::GetInstance();
		fecore.RegisterFactory(this);
	}
	FECoreBase* Create(FEModel* pfem) const { return new T(pfem); }
};

//-----------------------------------------------------------------------------
// Register a factory class
#define REGISTER_FECORE_FACTORY(theFactory) \
	static theFactory _##theFactory##_rc;

//-----------------------------------------------------------------------------
// Register a class using default creation parameters
#define REGISTER_FECORE_CLASS(theClass, ...) \
	static FERegisterClass_T<theClass> _##theClass##_rc(theClass::superClassID(), #theClass, theClass::BaseClassName(), __VA_ARGS__);

//-----------------------------------------------------------------------------
// Register a class using default creation parameters
#define REGISTER_FECORE_CLASS_EXPLICIT(theClass, theID, ...) \
	static FERegisterClass_T<theClass> _##theClass##_rc(theID, #theClass, theClass::BaseClassName(), __VA_ARGS__);

//-----------------------------------------------------------------------------
// version for classes that require template arguments
#define REGISTER_FECORE_CLASS_T(theClass, theArg, theName) \
	static FERegisterClass_T<theClass<theArg> > _##theClass##theArg##_rc(theClass<theArg>::superClassID(), #theClass, theClass<theArg>::BaseClassName(), theName);

//-----------------------------------------------------------------------------
// version for classes that require template arguments
#define REGISTER_FECORE_CLASS_T2(theClass, theArg1, theArg2, theName) \
	static FERegisterClass_T<theClass<theArg1, theArg2> > _##theClass##theArg1##theArg2##_rc(theClass<theArg1, theArg2>::superClassID(), #theClass, theClass<theArg1, theArg2>::BaseClassName(), theName);

//-----------------------------------------------------------------------------
// Create an instance of a class.
// This assumes that TBase is derived from FECoreBase and defines a class ID. 
template <typename TBase> inline TBase* fecore_new(const char* sztype, FEModel* pfem)
{
	FECoreKernel& fecore = FECoreKernel::GetInstance();
	return static_cast<TBase*>(fecore.Create(TBase::superClassID(), sztype, pfem));
//	return static_cast<TBase*>(fecore.Create(TBase::BaseClassName(), sztype, pfem));
}

//-----------------------------------------------------------------------------
// Create an instance of a class.
// This assumes that TBase is derived from FECoreBase and defines a class ID. 
template <typename TBase> inline TBase* fecore_new(int classIndex, FEModel* pfem)
{
	FECoreKernel& fecore = FECoreKernel::GetInstance();
	const FECoreFactory* f = fecore.GetFactoryClass(TBase::superClassID(), classIndex);
	if (f) return static_cast<TBase*>(f->Create(pfem));
	else return nullptr;
}

//-----------------------------------------------------------------------------
template <typename TClass> inline TClass* fecore_new_class(const char* szclass, FEModel* fem)
{
	int superId = TClass::superClassID();
	FECoreKernel& fecore = FECoreKernel::GetInstance();
	return static_cast<TClass*>(fecore.CreateClass(szclass, fem));
}

//-----------------------------------------------------------------------------
#define fecore_alloc(theClass, fem) fecore_new_class<theClass>(#theClass, fem)

//-----------------------------------------------------------------------------
// Three-parameter form of the fecore_new function for situations where the base class does not 
// define the classID value.
template <typename TBase> inline TBase* fecore_new(int sid, const char* sztype, FEModel* pfem)
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
	FEPluginFactory_T(const char* sz) : FECoreFactory(sid, nullptr, sz, nullptr){}
	FECoreBase* Create(FEModel* pfem) const { return new T(pfem); }
};

//------------------------------------------------------------------------------
// This is for functions exported from a plugin
#ifdef WIN32
#define FECORE_EXPORT extern "C" __declspec(dllexport)
#else
#define FECORE_EXPORT extern "C"
#endif
