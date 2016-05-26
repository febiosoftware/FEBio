#pragma once
#include <vector>
#include <assert.h>
#include "vec3d.h"
#include "vec2d.h"

//-----------------------------------------------------------------------------
class DumpStream;

//-----------------------------------------------------------------------------
enum FEDataType {
	FE_DOUBLE,
	FE_VEC2D,
	FE_VEC3D
};

//-----------------------------------------------------------------------------
class FEDataArray
{
public:
	//! default constructor
	FEDataArray(int dataType);

	//! copy constructor
	FEDataArray(const FEDataArray& map);

	//! assignment operator
	FEDataArray& operator = (const FEDataArray& map);

	//! allocate data
	bool resize(int nsize);

	//! return the size of array (units of data type)
	int size() const;

	//! get the value for a given facet index
	template <class T> T get(int n) const;

	//! set the value
	template <class T> bool set(const T& v);

	//! set the value
	template <class T> bool set(int n, const T& v);

	//! add a value
	template <class T> void push_back(const T& v);

public:
	//! get data type
	int DataType() const { return m_dataType; }

	//! return the buffer size (actual number of doubles)
	int BufferSize() const { return (int) m_val.size(); }

	//! serialization
	void Serialize(DumpStream& ar);

protected:
	int DataSize() const;

private:
	int	m_dataType;				//!< data type

	// default values
	double	m_defDouble;
	vec3d	m_defVec3d;
	vec2d	m_defVec2d;

	std::vector<double>	m_val;	//!< data values
};

template <> inline double FEDataArray::get<double>(int n) const
{
	assert(m_dataType == FE_DOUBLE);
	return	m_val[n];
}

template <> inline vec2d FEDataArray::get<vec2d>(int n) const
{
	assert(m_dataType == FE_VEC2D);
	return	vec2d(m_val[2*n], m_val[2*n+1]);
}

template <> inline vec3d FEDataArray::get<vec3d>(int n) const
{
	assert(m_dataType == FE_VEC3D);
	return	vec3d(m_val[3*n], m_val[3*n + 1], m_val[3*n+2]);
}

template <> inline void FEDataArray::push_back<double>(const double& v)
{
	assert(m_dataType == FE_DOUBLE);
	m_val.push_back(v);
}

template <> inline void FEDataArray::push_back<vec2d>(const vec2d& v)
{
	assert(m_dataType == FE_VEC2D);
	m_val.push_back(v.x());
	m_val.push_back(v.y());
}

template <> inline void FEDataArray::push_back<vec3d>(const vec3d& v)
{
	assert(m_dataType == FE_VEC3D);
	m_val.push_back(v.x);
	m_val.push_back(v.y);
	m_val.push_back(v.z);
}
