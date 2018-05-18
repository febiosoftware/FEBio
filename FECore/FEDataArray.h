#pragma once
#include <vector>
#include <assert.h>
#include <string>
#include "vec3d.h"
#include "vec2d.h"
#include "fecore_api.h"

//-----------------------------------------------------------------------------
class DumpStream;

//-----------------------------------------------------------------------------
enum FEDataType {
	FE_DOUBLE = 1,
	FE_VEC2D  = 2,
	FE_VEC3D  = 3
};

//-----------------------------------------------------------------------------
class FECORE_API FEDataArray
{
public:
	//! default constructor
	FEDataArray(int dataSize);

	virtual ~FEDataArray();

public:
	virtual void setValue(int n, double v) = 0;
	virtual void setValue(int n, const vec2d& v) = 0;
	virtual void setValue(int n, const vec3d& v) = 0;

	virtual void fillValue(double v) = 0;
	virtual void fillValue(const vec2d& v) = 0;
	virtual void fillValue(const vec3d& v) = 0;

protected:
	//! get the value for a given facet index
	template <class T> T get(int n) const;

	//! set the value
	template <class T> bool set(const T& v);

	//! set the value
	template <class T> bool set(int n, const T& v);

	//! add a value
	template <class T> void push_back(const T& v);

	//! allocate data
	bool resize(int nsize, double val = 0.0);

	//! set the data sized
	void SetDataSize(int dataSize);

public:
	//! get data size
	int DataSize() const { return m_dataSize; }

	// number of data items
	int DataCount() const { return m_dataCount; }

	//! return the buffer size (actual number of doubles)
	int BufferSize() const { return (int) m_val.size(); }

	//! serialization
	void Serialize(DumpStream& ar);

protected:
	//! copy constructor
	FEDataArray(const FEDataArray& map);

	//! assignment operator
	FEDataArray& operator = (const FEDataArray& map);

private:
	int	m_dataSize;			//!< size of each data item
	int	m_dataCount;		//!< number of data items

	std::vector<double>	m_val;	//!< data values
};

template <> inline double FEDataArray::get<double>(int n) const
{
	assert(m_dataSize == FE_DOUBLE);
	return	m_val[n];
}

template <> inline vec2d FEDataArray::get<vec2d>(int n) const
{
	assert(m_dataSize == FE_VEC2D);
	return	vec2d(m_val[2*n], m_val[2*n+1]);
}

template <> inline vec3d FEDataArray::get<vec3d>(int n) const
{
	assert(m_dataSize == FE_VEC3D);
	return	vec3d(m_val[3*n], m_val[3*n + 1], m_val[3*n+2]);
}

template <> inline void FEDataArray::push_back<double>(const double& v)
{
	assert(m_dataSize == FE_DOUBLE);
	m_val.push_back(v);
}

template <> inline void FEDataArray::push_back<vec2d>(const vec2d& v)
{
	assert(m_dataSize == FE_VEC2D);
	m_val.push_back(v.x());
	m_val.push_back(v.y());
}

template <> inline void FEDataArray::push_back<vec3d>(const vec3d& v)
{
	assert(m_dataSize == FE_VEC3D);
	m_val.push_back(v.x);
	m_val.push_back(v.y);
	m_val.push_back(v.z);
}

//-----------------------------------------------------------------------------
template <> inline bool FEDataArray::set<double>(int n, const double& v)
{
	assert(m_dataSize == FE_DOUBLE);
	m_val[n] = v;
	return true;
}

//-----------------------------------------------------------------------------
template <> inline bool FEDataArray::set<vec2d>(int n, const vec2d& v)
{
	assert(m_dataSize == FE_VEC2D);
	m_val[2 * n] = v.x();
	m_val[2 * n + 1] = v.y();
	return true;
}

//-----------------------------------------------------------------------------
template <> inline bool FEDataArray::set<vec3d>(int n, const vec3d& v)
{
	assert(m_dataSize == FE_VEC3D);
	m_val[3 * n] = v.x;
	m_val[3 * n + 1] = v.y;
	m_val[3 * n + 2] = v.z;
	return true;
}

//-----------------------------------------------------------------------------
template <> inline bool FEDataArray::set<double>(const double& v)
{
	assert(m_dataSize == FE_DOUBLE);
	for (int i = 0; i<(int)m_val.size(); ++i) m_val[i] = v;
	return true;
}

//-----------------------------------------------------------------------------
template <> inline bool FEDataArray::set<vec2d>(const vec2d& v)
{
	assert(m_dataSize == FE_VEC2D);
	for (int i = 0; i<(int)m_val.size(); i += 2)
	{
		m_val[i] = v.x();
		m_val[i + 1] = v.y();
	}
	return true;
}

//-----------------------------------------------------------------------------
template <> inline bool FEDataArray::set<vec3d>(const vec3d& v)
{
	assert(m_dataSize == FE_VEC3D);
	for (int i = 0; i<(int)m_val.size(); i += 3)
	{
		m_val[i] = v.x;
		m_val[i + 1] = v.y;
		m_val[i + 2] = v.z;
	}
	return true;
}
