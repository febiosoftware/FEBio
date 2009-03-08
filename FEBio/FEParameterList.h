#pragma once

#include <list>
#include <memory>
using namespace std;


//-----------------------------------------------------------------------------
// Different supported parameter types
enum FEParamType {
	FE_PARAM_INT,
	FE_PARAM_BOOL,
	FE_PARAM_DOUBLE,
	FE_PARAM_INTV = 100,
	FE_PARAM_DOUBLEV,
};

//-----------------------------------------------------------------------------
//! A material parameter
class FEParam
{
public:
	void*		m_pv;		// pointer to variable data
	FEParamType	m_itype;	// type of variable
	const char*	m_szname;	// name of the parameter
	int			m_ndim;		// dimension of array
	int			m_nlc;		// load curve number for dynamic parameters (-1 for static)

	FEParam()
	{
		m_pv = 0;
		m_itype = FE_PARAM_DOUBLE;
		m_ndim = 1;
		m_nlc = -1;
	}

public:
	//! retrieves the value for a non-array item
	template <class T> T& value() { return *((T*) m_pv); }

	//! retrieves the value for an array item
	template <class T> T* pvalue() { return (T*) m_pv; }
};

//-----------------------------------------------------------------------------
//! A list of material parameters
class FEParameterList
{
public:
	//! Add a parameter to the list
	void AddParameter(void* pv, FEParamType itype, int ndim, const char* sz);

	//! find a parameter using it's name
	FEParam* Find(const char* sz);

	//! returs the first parameter
	list<FEParam>::iterator first() { return m_pl.begin(); }

	//! number of parameters
	int Parameters() { return m_pl.size(); }

protected:
	list<FEParam>	m_pl;	//!< the actual parameter list
};

//-----------------------------------------------------------------------------
// To add parameter list to a material, simply do the following two steps
// 1) add the DECLARE_PARAMETER_LIST macro in the material class declaration
// 2) use the BEGIN_PARAMETER_LIST, ADD_PARAM and END_PARAMETER_LIST to
//    define a parameter list

// the following macro declares the parameter list for a material
#define DECLARE_PARAMETER_LIST() \
public: \
	virtual FEParameterList* GetParameterList();

// the BEGIN_PARAMETER_LIST defines the beginning of a parameter list
#define BEGIN_PARAMETER_LIST(theClass, baseClass) \
	FEParameterList* theClass::GetParameterList() { \
		FEParameterList* pl = baseClass::GetParameterList(); \

// the ADD_PARAMETER macro adds a parameter to the parameter list
#define ADD_PARAMETER(theParam, theType, theName) \
	pl->AddParameter(&theParam, theType, 1, theName);

// the ADD_PARAMETERV macro adds a parameter to the paramter list
// that is an array
#define ADD_PARAMETERV(theParam, theType, theDim, theName) \
	pl->AddParameter(theParam, theType, theDim, theName);

// the END_PARAMETER_LIST defines the end of a parameter list
#define END_PARAMETER_LIST() \
	return pl; }
