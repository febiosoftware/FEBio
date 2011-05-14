#pragma once
#include "FEBioFactory.h"
#include <vector>

class FEModel;

//-----------------------------------------------------------------------------
// This is the FEBio kernel class which manages the interactions between the 
// different executable modules. In particular, it manages the factory classes
// which are responsible for the life of different classes
class FEBioKernel
{
public:
	static FEBioKernel& GetInstance();

public:
	void RegisterClass(FEBioFactory* ptf) { m_Fac.push_back(ptf); }
	template <typename T> T* Create(const char* sztag, FEModel* pfem);

	template <typename T> const char* GetTypeStr(T* po);

protected:
	std::vector<FEBioFactory*>	m_Fac;

private:
	FEBioKernel(){}
	FEBioKernel(const FEBioKernel&){}
	static FEBioKernel* m_pKernel;
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
//! This class helps with the registration of a class with the FEBio framework
template <typename TDerived, typename TBase> class FERegisterClass_T : public FEBioFactory_T<TBase>
{
public:
	FERegisterClass_T(const char* sz) : FEBioFactory_T<TBase>(sz)
	{
		FEBioKernel& febio = FEBioKernel::GetInstance();
		febio.RegisterClass(this);
	}

	TDerived* Create(FEModel* pfem) { return new TDerived(pfem); }
	bool IsType(TBase* po) { return (dynamic_cast<TDerived*>(po) != 0); }
};

//-----------------------------------------------------------------------------
#define REGISTER_FEBIO_CLASS(theClass, theBase, theName) \
	static FERegisterClass_T<theClass, theBase> _##theClass##_rc(theName);
