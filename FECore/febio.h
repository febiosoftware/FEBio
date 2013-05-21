#pragma once
#include "FEBioFactory.h"
#include <vector>

//-----------------------------------------------------------------------------
// Forward declaration of model class
class FEModel;

//-----------------------------------------------------------------------------
// forward declaration of the log file
class Logfile;

//-----------------------------------------------------------------------------
// This is the FEBio kernel class that manages the interactions between the 
// different modules. In particular, it manages the factory classes
// which are responsible for the creation of different classes that are registered
// with the kernel.

// TODO: I am using template definitions in this class which means that different plugins
//       will see a different interface to FEBioKernel, depending on the classes that
//		 are registered in the pluing. Could this be a problem?
//
class FEBioKernel
{
public:
	// Do not call this function from a plugin as it will not return the correct
	// instance. Instead, use the FEBioKernel object that is passed in the RegisterFEBioPlugin method
	static FEBioKernel& GetInstance();

	// Get the logfile
	Logfile& GetLogfile();

public:
	void RegisterClass(FEBioFactory* ptf) { m_Fac.push_back(ptf); }

public:
	template <typename T> T* Create(const char* sztag, FEModel* pfem);

	template <typename T> const char* GetTypeStr(T* po);

	template <typename T> const char* GetTypeStr(int i);

	template <typename T> int Count();

private:
	std::vector<FEBioFactory*>	m_Fac;	// list of registered factory classes

	Logfile*	m_plog;	// keep a pointer to the logfile (used by plugins)

private: // make singleton
	FEBioKernel();
	FEBioKernel(const FEBioKernel&){}
	void operator = (const FEBioKernel&){}

private:
	static FEBioKernel* m_pKernel;	// the one-and-only kernel object
};

//-----------------------------------------------------------------------------
template <typename T> inline T* FEBioKernel::Create(const char* sztag, FEModel* pfem)
{
	std::vector<FEBioFactory*>::iterator pf;
	for (pf=m_Fac.begin(); pf!= m_Fac.end(); ++pf)
	{
		FEBioFactory_T<T>* pfac = dynamic_cast<FEBioFactory_T<T>*>(*pf);
		if (pfac)
		{
			if (strcmp(pfac->GetTypeStr(), sztag) == 0) return pfac->Create(pfem);
		}
	}
	return 0;
}

//-----------------------------------------------------------------------------
template <typename T> inline const char* FEBioKernel::GetTypeStr(T* po)
{
	std::vector<FEBioFactory*>::iterator pf;
	for (pf=m_Fac.begin(); pf!= m_Fac.end(); ++pf)
	{
		FEBioFactory_T<T>* pfac = dynamic_cast<FEBioFactory_T<T>*>(*pf);
		if (pfac)
		{
			if (pfac->IsType(po) == true) return pfac->GetTypeStr();
		}
	}
	return 0;
}

//-----------------------------------------------------------------------------
template <typename T> inline const char* FEBioKernel::GetTypeStr(int n)
{
	int i = 0;
	std::vector<FEBioFactory*>::iterator pf;
	for (pf=m_Fac.begin(); pf!= m_Fac.end(); ++pf)
	{
		FEBioFactory_T<T>* pfac = dynamic_cast<FEBioFactory_T<T>*>(*pf);
		if (pfac)
		{
			if (i == n) return pfac->GetTypeStr();
			++i;
		}
	}
	return 0;
}

//-----------------------------------------------------------------------------
template <typename T> inline int FEBioKernel::Count()
{
	int N = 0;
	std::vector<FEBioFactory*>::iterator pf;
	for (pf=m_Fac.begin(); pf!= m_Fac.end(); ++pf)
	{
		FEBioFactory_T<T>* pfac = dynamic_cast<FEBioFactory_T<T>*>(*pf);
		if (pfac) N++;
	}
	return N;
}

//-----------------------------------------------------------------------------
//! This class helps with the registration of a class with the FEBio framework
//! The TDerived class is the typename of the class to be registered.
//! The TBase class is the name of the base class, that is the class that TDerived inherits from
template <typename TDerived, typename TBase> class FERegisterClass_T : public FEBioFactory_T<TBase>
{
public:
	FERegisterClass_T(const char* sz) : FEBioFactory_T<TBase>(sz)
	{
		FEBioKernel& febio = FEBioKernel::GetInstance();
		febio.RegisterClass(this);
	}

	TBase* Create(FEModel* pfem) { return new TDerived(pfem); }
	bool IsType(TBase* po) { return (dynamic_cast<TDerived*>(po) != 0); }
};

//-----------------------------------------------------------------------------
// Register a class using default creation parameters
#define REGISTER_FEBIO_CLASS(theClass, theBase, theName) \
	static FERegisterClass_T<theClass, theBase> _##theClass##_rc(theName);

//-----------------------------------------------------------------------------
// version for classes that require template arguments
#define REGISTER_FEBIO_CLASS_T(theClass, theBase, theArg, theName) \
	static FERegisterClass_T<theClass<theArg>, theBase> _##theClass##theArg##_rc(theName);
