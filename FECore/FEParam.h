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
	FE_PARAM_FUNC1D,
	FE_PARAM_TENS3DRS,
	FE_PARAM_STD_VECTOR_DOUBLE,
	FE_PARAM_DOUBLE_MAPPED,
	FE_PARAM_VEC3D_MAPPED
};

//-----------------------------------------------------------------------------
// class describing the value of parameter
class FECORE_API FEParamValue
{
private:
	void*			m_pv;		// pointer to variable data
	FEParamType		m_itype;	// type of variable
	int				m_ndim;		// dimension of array (or 1 for non-array types)

public:

	FEParamValue()
	{
		m_pv = 0;
		m_itype = FE_PARAM_INVALID;
		m_ndim = -1;
	}

	explicit FEParamValue(double& v)
	{
		m_pv = &v;
		m_itype = FE_PARAM_DOUBLE;
		m_ndim = 1;
	}

	explicit FEParamValue(vec3d& v)
	{
		m_pv = &v;
		m_itype = FE_PARAM_VEC3D;
		m_ndim = 1;
	}

	explicit FEParamValue(void* data, FEParamType itype, int dim = 1)
	{
		m_pv = data;
		m_itype = itype;
		m_ndim = dim;
	}

	bool isValid() const { return (m_pv != 0); }

	FEParamType type() const { return m_itype; }

	void* data_ptr() const { return m_pv; }

	int dim() const
	{
		return m_ndim;
	}

	template <typename T> T& value() { return *((T*)m_pv); }
	template <typename T> const T& value() const { return *((T*)m_pv); }

	void Serialize(DumpStream& ar);
};

//-----------------------------------------------------------------------------
//! This class describes a user-defined parameter
class FECORE_API FEParam
{
private:
	FEParamValue	m_val;	// stores the value of the parameter

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

	// get the value
	FEParamValue& paramValue() { return m_val; }

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
	int dim() const { return m_val.dim(); }

	// parameter type
	FEParamType type() const { return m_val.type(); }

	// data pointer
	void* data_ptr() const { return m_val.data_ptr(); }

	// set the load curve ID and scale factor
	void SetLoadCurve(int lc);
	void SetLoadCurve(int lc, double s);
	void SetLoadCurve(int lc, const vec3d& v );

	// get the load curve ID (or -1 if none)
	int GetLoadCurve() const { return m_nlc; }

	// get the scale factors
	FEParamValue GetScale() 
	{ 
		if (m_val.type() == FE_PARAM_DOUBLE)
			return FEParamValue(m_scl); 
		else if (m_val.type() == FE_PARAM_VEC3D)
			return FEParamValue(m_vscl);
		else return FEParamValue();
	}
	double& GetScaleDouble() { return m_scl; }
	vec3d& GetScaleVec3d () { return m_vscl; }

	// Copy the state of one parameter to this parameter.
	// This requires that the parameters are compatible (i.e. same type, etc.)
	bool CopyState(const FEParam& p);

	void setParent(FEParamContainer* pc) { m_parent = pc; }
	FEParamContainer* parent() { return m_parent; }

public:
	void Serialize(DumpStream& ar);

public:
	//! retrieves the value for a non-array item
	template <class T> T& value() { return *((T*)m_val.data_ptr()); }

	//! retrieves the value for a non-array item
	template <class T> const T& value() const { return *((T*)m_val.data_ptr()); }

	//! retrieves the value for an array item
	template <class T> T* pvalue() { return (T*)m_val.data_ptr(); }

	//! retrieves the value for an array item
	template <class T> T value(int i) const { return ((T*)m_val.data_ptr())[i]; }

	//! retrieves pointer to element in array
	template <class T> T* pvalue(int n);

	//! override the template for char pointers
	char* cvalue() { return (char*)m_val.data_ptr(); }
};

//-----------------------------------------------------------------------------
//! Retrieves a pointer to element in array
template<class T> inline T* FEParam::pvalue(int n)
{
	assert((n >= 0) && (n < m_val.dim()));
	return &(pvalue<T>()[n]);
}

//-----------------------------------------------------------------------------
FEParamValue GetParameterComponent(const ParamString& paramName, FEParam* param);
