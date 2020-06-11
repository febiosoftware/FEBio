/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2020 University of Utah, The Trustees of Columbia University in 
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
#include <string.h>
#include <stdio.h>
#include "version.h"

//-----------------------------------------------------------------------------
// Forward declarations
class FEModel;
class LinearSolver;

//-----------------------------------------------------------------------------
//! This is the FECore kernel class that manages the interactions between the 
//! different modules. In particular, it manages the factory classes
//! which are responsible for the creation of different classes that are registered
//! with the kernel.
class FECORE_API FECoreKernel
{
	enum { ALL_MODULES = 0xFFFF };

	struct Module
	{
		const char*		szname;	// name of module
		unsigned int	id;		// unqiue ID (this is a bit value)
		unsigned int	flags;	// ID + IDs of dependent modules
		int				m_alloc_id;	// ID of allocator
	};

public:
	// Do not call this function from a plugin as it will not return the correct
	// instance. Instead, use the FECoreKernel object that is passed in the PluginInitialize method
	static FECoreKernel& GetInstance();

	// set the instance of the kernel
	static void SetInstance(FECoreKernel* pkernel);

public:
	//! Register a class with the framework
	void RegisterFactory(FECoreFactory* ptf);

	//! Create a specific using a superclass ID and an alias
	void* Create(int superClassID, const char* szalias, FEModel* pfem);

	//! Create a specific class
	void* CreateClass(const char* szclassName, FEModel* fem);

	//! Create a class from a class descriptor
	void* Create(int superClassID, FEModel* pfem, const ClassDescriptor& cd);

	//! count the number of registered classes with a given super-class id
	int Count(SUPER_CLASS_ID sid);

	//! List the registered classes with a given super-class id
	void List(SUPER_CLASS_ID sid);

	//! Get the number of registered factory classes
	int FactoryClasses();

	//! return a factory class
	const FECoreFactory* GetFactoryClass(int i);

	//! return a factory class
	const FECoreFactory* GetFactoryClass(int classID, int i);

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

public: // Modules

	//! count modules
	int Modules() const;

	//! Create a module (also makes it the active module)
	bool CreateModule(const char* szmodule);

	//! set the active module
	bool SetActiveModule(const char* szmodule);

	//! set a dependency on a module
	bool SetModuleDependency(const char* szmodule);

	//! remove a module
	bool RemoveModule(const char* szmodule);

	//! Get a module
	const char* GetModuleName(int i) const;

	//! set the spec ID. Features with a matching spec ID will be preferred
	//! set spec ID to -1 to stop caring
	void SetSpecID(int nspec);

public:
	//! Register a new domain class
	void RegisterDomain(FEDomainFactory* pf);

	//! Create a domain of a certain type
	FEDomain* CreateDomain(const FE_Element_Spec& spec, FEMesh* pm, FEMaterial* pmat);

public:
	//! set the default linear solver
	FECoreFactory* SetDefaultSolverType(const char* sztype);

	void SetDefaultSolver(ClassDescriptor* linsolve);

	//! get the linear solver type
	const char* GetLinearSolverType() const;

	//! Get a linear solver
	LinearSolver* CreateDefaultLinearSolver(FEModel* fem);

private:
	std::vector<FECoreFactory*>			m_Fac;	// list of registered factory classes
	std::vector<FEDomainFactory*>		m_Dom;	// list of domain factory classes

	std::string			m_default_solver_type;	// default linear solver
	ClassDescriptor*	m_default_solver;

	// module list
	vector<Module>	m_modules;
	int				m_activeModule;

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
	FERegisterClass_T(SUPER_CLASS_ID sid, const char* szclass, const char* szalias, int spec = -1) : FECoreFactory(sid, szclass, szalias, spec)
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
	static FERegisterClass_T<theClass> _##theClass##_rc(theClass::classID(), #theClass, __VA_ARGS__);

//-----------------------------------------------------------------------------
// Register a class using default creation parameters
#define REGISTER_FECORE_CLASS_EXPLICIT(theClass, theID, ...) \
	static FERegisterClass_T<theClass> _##theClass##_rc(theID, #theClass, __VA_ARGS__);

//-----------------------------------------------------------------------------
// version for classes that require template arguments
#define REGISTER_FECORE_CLASS_T(theClass, theArg, theName) \
	static FERegisterClass_T<theClass<theArg> > _##theClass##theArg##_rc(theClass<theArg>::classID(), 0, theName);

//-----------------------------------------------------------------------------
// Create an instance of a class.
// This assumes that TBase is derived from FECoreBase and defines a class ID. 
template <typename TBase> inline TBase* fecore_new(const char* sztype, FEModel* pfem)
{
	FECoreKernel& fecore = FECoreKernel::GetInstance();
	return static_cast<TBase*>(fecore.Create(TBase::classID(), sztype, pfem));
}

//-----------------------------------------------------------------------------
// Create an instance of a class.
// This assumes that TBase is derived from FECoreBase and defines a class ID. 
template <typename TBase> inline TBase* fecore_new(int classIndex, FEModel* pfem)
{
	FECoreKernel& fecore = FECoreKernel::GetInstance();
	const FECoreFactory* f = fecore.GetFactoryClass(TBase::classID(), classIndex);
	if (f) return static_cast<TBase*>(f->Create(pfem));
	else return nullptr;
}

//-----------------------------------------------------------------------------
template <typename TClass> inline TClass* fecore_new_class(const char* szclass, FEModel* fem)
{
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
	FEPluginFactory_T(const char* sz) : FECoreFactory(sid, nullptr, sz){}
	FECoreBase* Create(FEModel* pfem) const { return new T(pfem); }
};

//------------------------------------------------------------------------------
// This is for functions exported from a plugin
#ifdef WIN32
#define FECORE_EXPORT extern "C" __declspec(dllexport)
#else
#define FECORE_EXPORT extern "C"
#endif
