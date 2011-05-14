#pragma once

//-----------------------------------------------------------------------------
class FEModel;

//-----------------------------------------------------------------------------
class FEBioFactory
{
public:
	virtual ~FEBioFactory(){}

	// return the type string identifier
	const char* GetTypeStr() { return m_sztype; }

protected:
	// constructor
	FEBioFactory(const char* sztype) { m_sztype = sztype; }

protected:
	const char*	m_sztype;	//!< class type string
};

//-----------------------------------------------------------------------------
// The FEBioFactory_T class is a template class that can be used to create 
// factory classes for FEBio. Factory are used by the framework to create instance
// of FEBio classes. 
template <typename T> class FEBioFactory_T : public FEBioFactory
{
protected:
	// contructor for the factory
	FEBioFactory_T(const char* sztype) : FEBioFactory(sztype){}

public:
	// create an instance of TDerived class
	virtual T* Create(FEModel* pfem) = 0;

	// check type of class
	virtual bool IsType(T* po) = 0;
};
