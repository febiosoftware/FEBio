#pragma once

#include "vec3d.h"
#include <assert.h>
#include <list>
#include <memory>
#include "DumpFile.h"
using namespace std;


//-----------------------------------------------------------------------------
// Different supported parameter types
enum FEParamType {
	FE_PARAM_INT,
	FE_PARAM_BOOL,
	FE_PARAM_DOUBLE,
	FE_PARAM_VEC3D,
	FE_PARAM_MAT3D,
	FE_PARAM_MAT3DS,
	FE_PARAM_IMAGE_3D,
	FE_PARAM_STRING,
	FE_PARAM_INTV = 100,
	FE_PARAM_DOUBLEV,
};

//-----------------------------------------------------------------------------
// Range for parameters
enum FEParamRange {
	FE_DONT_CARE,			// parameter can take on any value
	FE_GREATER,				// parameter must be greater than val (> val)
	FE_GREATER_OR_EQUAL,	// parameter must be greater or equal than val ( >= val)
	FE_LESS,				// parameter must be less than val ( < val)
	FE_LESS_OR_EQUAL,		// parameter must be less or equal than val ( <= val)
	FE_OPEN,				// parameter must be inside the open interval (min,max)
	FE_CLOSED,				// parameter must be inside the closed interval [min, max]
	FE_LEFT_OPEN,			// parameter must be inside the half-open interval (min, max]
	FE_RIGHT_OPEN,			// parameter must be inside the right-open interval [min, max)
	FE_NOT_EQUAL,			// parameter must not equal val ( != val)
};

//-----------------------------------------------------------------------------
//! This class describes a user-defined parameter
class FEParam
{
public:
	void*		m_pv;		// pointer to variable data
	FEParamType	m_itype;	// type of variable
	const char*	m_szname;	// name of the parameter
	int			m_ndim;		// dimension of array
	int			m_nlc;		// load curve number for dynamic parameters (-1 for static)
	double		m_scl;		// load curve scale factor
	vec3d		m_vscl;		// scale factor for vectors

	// range stuff
	FEParamRange	m_irange;
	union
	{
		// this value is used as the minimum for a range parameter
		// or for the inf/sup for a half-open or half-closed parameter
		int		m_imin;
		double	m_dmin;
	};
	union
	{
		// this value is used as the maximum for a range parameter
		int		m_imax;
		double	m_dmax;
	};

	FEParam()
	{
		m_pv = 0;
		m_itype = FE_PARAM_DOUBLE;
		m_irange = FE_DONT_CARE;
		m_ndim = 1;
		m_nlc = -1;
		m_scl = 1.0;
		m_vscl = vec3d(0,0,0);
	}

	bool is_inside_range();

public:
	//! retrieves the value for a non-array item
	template <class T> T& value() { return *((T*) m_pv); }

	//! retrieves the value for an array item
	template <class T> T* pvalue() { return (T*) m_pv; }

	//! retrieves pointer to element in array
	template <class T> T* pvalue(int n);

	//! assignment operators
	void operator = (double g) { assert(m_itype == FE_PARAM_DOUBLE); value<double>() = g; }
	void operator = (int    n) { assert(m_itype == FE_PARAM_INT   ); value<int   >() = n; }
	void operator = (bool   b) { assert(m_itype == FE_PARAM_BOOL  ); value<bool  >() = b; }
	void operator = (vec3d  v) { assert(m_itype == FE_PARAM_VEC3D ); value<vec3d >() = v; }
	void operator = (mat3d  m) { assert(m_itype == FE_PARAM_MAT3D ); value<mat3d >() = m; }
	void operator = (mat3ds m) { assert(m_itype == FE_PARAM_MAT3DS); value<mat3ds>() = m; }

	//! override the template for char pointers
	char* cvalue() { return (char*) m_pv; }
};

//-----------------------------------------------------------------------------
//! Retrieves a pointer to element in array
//! \todo only works with doubles for now
template<class T> inline T* FEParam::pvalue(int n)
{
	assert((m_itype==FE_PARAM_DOUBLE)||(m_itype==FE_PARAM_DOUBLEV));
	assert((n >= 0) && (n < m_ndim));
	return &(pvalue<double>()[n]);
}

//-----------------------------------------------------------------------------
// helper class for defining ranges
class RANGE
{
public:
	FEParamRange	m_rt;
	double			m_fmin;
	double			m_fmax;
public:
	RANGE(FEParamRange rt) { m_rt = rt; m_fmin = m_fmax = 0.0; }
	RANGE(FEParamRange rt, double fval) { m_rt = rt; m_fmin = fval; m_fmax = 0.0; }
	RANGE(FEParamRange rt, double fmin, double fmax) { m_rt = rt; m_fmin = fmin; m_fmax = fmax; }
};

inline RANGE FE_RANGE_DONT_CARE() { return RANGE(FE_DONT_CARE); }
inline RANGE FE_RANGE_GREATER         (double f) { return RANGE(FE_GREATER         , f); }
inline RANGE FE_RANGE_GREATER_OR_EQUAL(double f) { return RANGE(FE_GREATER_OR_EQUAL, f); }
inline RANGE FE_RANGE_LESS            (double f) { return RANGE(FE_LESS            , f); }
inline RANGE FE_RANGE_LESS_OR_EQUAL   (double f) { return RANGE(FE_LESS_OR_EQUAL   , f); }
inline RANGE FE_RANGE_OPEN      (double fmin, double fmax) {return RANGE(FE_OPEN      , fmin, fmax); }
inline RANGE FE_RANGE_CLOSED    (double fmin, double fmax) {return RANGE(FE_CLOSED    , fmin, fmax); }
inline RANGE FE_RANGE_LEFT_OPEN (double fmin, double fmax) {return RANGE(FE_LEFT_OPEN , fmin, fmax); }
inline RANGE FE_RANGE_RIGHT_OPEN(double fmin, double fmax) {return RANGE(FE_RIGHT_OPEN, fmin, fmax); }
inline RANGE FE_RANGE_NOT_EQUAL (double f) { return RANGE(FE_NOT_EQUAL, f); }

//-----------------------------------------------------------------------------
// forward declaration of the param container
class FEParamContainer;

//-----------------------------------------------------------------------------
//! A list of material parameters
class FEParameterList
{
public:
	FEParameterList(FEParamContainer* pc) : m_pc(pc) {}
	virtual ~FEParameterList(){}

	//! assignment operator
	void operator = (FEParameterList& l);

	//! Add a parameter to the list
	void AddParameter(void* pv, FEParamType itype, int ndim, const char* sz);

	//! Add a parameter to the list
	void AddParameter(void* pv, FEParamType type, int ndim, FEParamRange rng, double fmin, double fmax, const char* sz);

	//! find a parameter using it's name (the safe way)
	FEParam* Find(const char* sz);

	//! get a parameter (the dangerous way)
	FEParam& operator [] (const char* sz) { return *Find(sz); }

	//! returs the first parameter
	list<FEParam>::iterator first() { return m_pl.begin(); }

	//! number of parameters
	int Parameters() { return m_pl.size(); }

	//! return the parameter container
	FEParamContainer* GetContainer() { return m_pc; }

protected:
	FEParamContainer*	m_pc;	//!< parent container
	list<FEParam>		m_pl;	//!< the actual parameter list
};

//-----------------------------------------------------------------------------
// helper class for retrieving parameters
class ParamString
{
public:
	//! constructor
	explicit ParamString(const char* sz);

	//! copy constructor
	ParamString(const ParamString& p);

	//! assignment operator
	void operator = (const ParamString& p);

	//! destructor
	~ParamString();

	//! see how many components the string has
	int count() const;

	//! return the next component
	ParamString next() const;

	//! return the next name in the list
	const char* c_str() const { return m_sz; }

	//! compare to a string
	bool operator == (const char* sz) const;

private:
	explicit ParamString(char* sz, int nc, int nl) : m_sz(sz), m_nc(nc), m_nl(nl) {}

private:
	char*	m_sz;	//!< string containing parameter names
	int		m_nc;	//!< number of string components
	int		m_nl;	//!< total string length
};

//-----------------------------------------------------------------------------
//! Base class for classes that wish to support parameter lists
class FEParamContainer
{
public:
	//! constructor
	FEParamContainer();

	//! destructor
	virtual ~FEParamContainer();

	//! return the material's parameter list
	FEParameterList& GetParameterList();

	//! find a parameter using it's name
	virtual FEParam* GetParameter(const ParamString& s);

	//! serialize parameter data
	virtual void Serialize(DumpFile& ar);

public:
	//! This function is called after the parameter was read in from the input file.
	//! It can be used to do additional processing when a parameter is read in.
	virtual void SetParameter(FEParam& p) {}

	//! If a parameter has attributes, this function will be called
	virtual bool SetParameterAttribute(FEParam& p, const char* szatt, const char* szval) { return false; }

protected:
	//! This function will be overridden by each class that defines a parameter list
	virtual void BuildParamList() {}

	//! Add a parameter to the list
	void AddParameter(void* pv, FEParamType itype, int ndim, const char* sz);

	//! Add a parameter to the list
	void AddParameter(void* pv, FEParamType type, int ndim, RANGE rng, const char* sz);

private:
	FEParameterList*	m_pParam;	//!< parameter list
};

//-----------------------------------------------------------------------------
// To add parameter list to a material, simply do the following two steps
// 1) add the DECLARE_PARAMETER_LIST macro in the material class declaration
// 2) use the BEGIN_PARAMETER_LIST, ADD_PARAM and END_PARAMETER_LIST to
//    define a parameter list

// the following macro declares the parameter list for a material
#define DECLARE_PARAMETER_LIST() \
protected: \
	virtual void BuildParamList();

// the BEGIN_PARAMETER_LIST defines the beginning of a parameter list
#define BEGIN_PARAMETER_LIST(theClass, baseClass) \
	void theClass::BuildParamList() { \
			baseClass::BuildParamList(); \

// the ADD_PARAMETER macro adds a parameter to the parameter list
#define ADD_PARAMETER(theParam, theType, theName) \
	AddParameter(&theParam, theType, 1, theName);

// the ADD_PARAMETERV macro adds a parameter to the paramter list
// that is an array
#define ADD_PARAMETERV(theParam, theType, theDim, theName) \
	AddParameter(theParam, theType, theDim, theName);

// the ADD_PARAMETER2 macro adds a parameter with range checking to the parameter list
#define ADD_PARAMETER2(theParam, theType, theRange, theName) \
	AddParameter(&theParam, theType, 1, theRange, theName);

// the END_PARAMETER_LIST defines the end of a parameter list
#define END_PARAMETER_LIST() \
	}
