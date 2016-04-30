#pragma once
#include "vec3d.h"
#include "mat3d.h"
#include <assert.h>

//-----------------------------------------------------------------------------
class FEParamValidator;

//-----------------------------------------------------------------------------
// Different supported parameter types
enum FEParamType {
	FE_PARAM_INT,
	FE_PARAM_BOOL,
	FE_PARAM_DOUBLE,
	FE_PARAM_VEC2D,
	FE_PARAM_VEC3D,
	FE_PARAM_MAT3D,
	FE_PARAM_MAT3DS,
	FE_PARAM_IMAGE_3D,
	FE_PARAM_STRING,
	FE_PARAM_SURFACE_MAP
};

//-----------------------------------------------------------------------------
//! This class describes a user-defined parameter
class FEParam
{
public:
	// TODO: I want to look into the idea of generalizing this to "controllers". 
	//       A controller would be anything that affects a parameter's value. Right now,
	//       only load curves are used, but other controllers can be a mathematical expression,
	//       or even a user-controlled interactive controller. 
	int			m_nlc;		// load curve number for dynamic parameters (-1 for static)

	// Can I put these two variables in a union?
	double		m_scl;		// load curve scale factor
	vec3d		m_vscl;		// scale factor for vectors

private:
	void*		m_pv;		// pointer to variable data
	FEParamType	m_itype;	// type of variable
	int			m_ndim;		// dimension of array
	const char*	m_szname;	// name of the parameter

	// parameter validator
	FEParamValidator*	m_pvalid;

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

	// parameter dimension
	int dim() const { return m_ndim; }

	// parameter type
	FEParamType type() const { return m_itype; }

	// data pointer
	void* data_ptr() const { return m_pv; }

public:
	//! retrieves the value for a non-array item
	template <class T> T& value() { return *((T*) m_pv); }

	//! retrieves the value for a non-array item
	template <class T> T value() const { return *((T*)m_pv); }

	//! retrieves the value for an array item
	template <class T> T* pvalue() { return (T*) m_pv; }

	//! retrieves the value for an array item
	template <class T> T value(int i) const { return ((T*)m_pv)[i]; }

	//! retrieves pointer to element in array
	template <class T> T* pvalue(int n);

	//! assignment operators
	// NOTE: This doesn't work if a load curve is used (use setvalue instead).
	void operator = (double g) { assert(m_itype == FE_PARAM_DOUBLE); value<double>() = g; }
	void operator = (int    n) { assert(m_itype == FE_PARAM_INT   ); value<int   >() = n; }
	void operator = (bool   b) { assert(m_itype == FE_PARAM_BOOL  ); value<bool  >() = b; }
	void operator = (vec2d  v) { assert(m_itype == FE_PARAM_VEC2D ); value<vec2d >() = v; }
	void operator = (vec3d  v) { assert(m_itype == FE_PARAM_VEC3D ); value<vec3d >() = v; }
	void operator = (mat3d  m) { assert(m_itype == FE_PARAM_MAT3D ); value<mat3d >() = m; }
	void operator = (mat3ds m) { assert(m_itype == FE_PARAM_MAT3DS); value<mat3ds>() = m; }

	// only implemented for double parameters
	void setvalue(double v)
	{
		if (m_nlc == -1) value<double>() = v;
		else
		{
			assert(m_itype == FE_PARAM_DOUBLE);
			m_scl = v;
		}
	}

	//! override the template for char pointers
	char* cvalue() { return (char*) m_pv; }
};

//-----------------------------------------------------------------------------
//! Retrieves a pointer to element in array
template<class T> inline T* FEParam::pvalue(int n)
{
	assert((n >= 0) && (n < m_ndim));
	return &(pvalue<T>()[n]);
}
