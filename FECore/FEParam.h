#pragma once
#include "vec3d.h"
#include "mat3d.h"
#include <assert.h>
#include <vector>
#include "fecore_api.h"
#include "ParamString.h"

//-----------------------------------------------------------------------------
class FEParamValidator;
class DumpStream;
class FEParamContainer;

//-----------------------------------------------------------------------------
// Different supported parameter types
enum FEParamType {
	FE_PARAM_INVALID,
	FE_PARAM_INT,
	FE_PARAM_BOOL,
	FE_PARAM_DOUBLE,
	FE_PARAM_VEC2D,
	FE_PARAM_VEC3D,
	FE_PARAM_MAT3D,
	FE_PARAM_MAT3DS,
	FE_PARAM_IMAGE_3D,
	FE_PARAM_STRING,
	FE_PARAM_DATA_ARRAY,
	FE_PARAM_TENS3DRS,
	FE_PARAM_STD_STRING,
	FE_PARAM_STD_VECTOR_DOUBLE,
	FE_PARAM_STD_VECTOR_VEC2D,
	FE_PARAM_DOUBLE_MAPPED,
	FE_PARAM_VEC3D_MAPPED
};

//-----------------------------------------------------------------------------
// Parameter flags
enum FEParamFlag {
	FE_PARAM_ATTRIBUTE = 0x01		// parameter will be read as attribute
};

class FEParam;

//-----------------------------------------------------------------------------
// class describing the value of parameter
class FECORE_API FEParamValue
{
private:
	void*			m_pv;		// pointer to variable data
	FEParamType		m_itype;	// type of variable (this is not the type of the param!)
	FEParam*		m_param;	// the parameter
	int				m_index;	// index of paramter value

public:

	FEParamValue()
	{
		m_pv = 0;
		m_itype = FE_PARAM_INVALID;
		m_param = 0;
		m_index = -1;
	}

	explicit FEParamValue(FEParam* p, void* v, FEParamType itype, int index = 0)
	{
		m_pv = v;
		m_itype = itype;
		m_param = p;
		m_index = index;
	}

	bool isValid() const { return (m_pv != 0); }

	FEParamType type() const { return m_itype; }

	void* data_ptr() const { return m_pv; }

	int index() const { return m_index; }

	FEParam* param() { return m_param; }

	template <typename T> T& value() { return *((T*)m_pv); }
	template <typename T> const T& value() const { return *((T*)m_pv); }
};

//-----------------------------------------------------------------------------
//! This class describes a user-defined parameter
class FECORE_API FEParam
{
private:
	void*			m_pv;		// pointer to variable data
	int				m_dim;		// dimension (in case data is array)
	FEParamType		m_type;		// type of variable
	unsigned int	m_flag;		// parameter flags

	const char*	m_szname;	// name of the parameter
	const char*	m_szenum;	// enumerate values for ints

	// parameter validator
	FEParamValidator*	m_pvalid;

	// TODO: I want to look into the idea of generalizing this to "controllers". 
	//       A controller would be anything that affects a parameter's value. Right now,
	//       only load curves are used, but other controllers can be a mathematical expression,
	//       or even a user-controlled interactive controller. 
	int			m_nlc;		// load curve number for dynamic parameters (-1 for static)

	// Can I put these two variables in a union?
	double		m_scl;		// load curve scale factor
	vec3d		m_vscl;		// scale factor for vectors

	FEParamContainer* m_parent;	// parent object of parameter

public:
	// constructor
	FEParam(void* pdata, FEParamType itype, int ndim, const char* szname);
	FEParam(const FEParam& p);
	FEParam& operator = (const FEParam& p);

	// set the parameter's validator
	void SetValidator(FEParamValidator* pvalid);

	// see if the parameter's value is valid
	bool is_valid() const;

	// return the name of the parameter
	const char* name() const { return m_szname; }

	// return the enum values
	const char* enums() const { return m_szenum; }

	// set the enum values (\0 separated. Make sure the end of the string has two \0's)
	void SetEnums(const char* sz) { m_szenum = sz; }

	// parameter dimension
	int dim() const { return m_dim; }

	// parameter type
	FEParamType type() const { return m_type; }

	// data pointer
	void* data_ptr() const { return m_pv; }

	// set the load curve ID and scale factor
	void SetLoadCurve(int lc);
	void SetLoadCurve(int lc, double s);
	void SetLoadCurve(int lc, const vec3d& v );

	// get the load curve ID (or -1 if none)
	int GetLoadCurve() const { return m_nlc; }

	// get the param value
	FEParamValue paramValue(int i = -1);

	// get the scale factors
	double& GetScaleDouble() { return m_scl; }
	vec3d& GetScaleVec3d () { return m_vscl; }

	// Copy the state of one parameter to this parameter.
	// This requires that the parameters are compatible (i.e. same type, etc.)
	bool CopyState(const FEParam& p);

	void setParent(FEParamContainer* pc) { m_parent = pc; }
	FEParamContainer* parent() { return m_parent; }

	void SetFlags(unsigned int flags) { m_flag = flags; }
	unsigned int GetFlags() const { return m_flag; }

public:
	void Serialize(DumpStream& ar);

public:
	//! retrieves the value for a non-array item
	template <class T> T& value() { return *((T*) data_ptr()); }

	//! retrieves the value for a non-array item
	template <class T> const T& value() const { return *((T*) data_ptr()); }

	//! retrieves the value for an array item
	template <class T> T* pvalue() { return (T*) data_ptr(); }

	//! retrieves the value for an array item
	template <class T> T& value(int i) { return ((T*)data_ptr())[i]; }
	template <class T> T value(int i) const { return ((T*) data_ptr())[i]; }

	//! retrieves pointer to element in array
	template <class T> T* pvalue(int n);

	//! override the template for char pointers
	char* cvalue() { return (char*) data_ptr(); }
};

//-----------------------------------------------------------------------------
//! Retrieves a pointer to element in array
template<class T> inline T* FEParam::pvalue(int n)
{
	assert((n >= 0) && (n < m_dim));
	return &(pvalue<T>()[n]);
}

//-----------------------------------------------------------------------------
FEParamValue GetParameterComponent(const ParamString& paramName, FEParam* param);
