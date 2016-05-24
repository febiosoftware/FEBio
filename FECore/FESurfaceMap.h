#pragma once
#include <vector>
#include <string>
#include <FECore/vec3d.h>
#include <assert.h>

//-----------------------------------------------------------------------------
class FESurface;
class FEFacetSet;
class DumpStream;

//-----------------------------------------------------------------------------
typedef int FEFacetIndex;

//-----------------------------------------------------------------------------
enum FEDataType {
	FE_DOUBLE,
	FE_VEC3D
};

//-----------------------------------------------------------------------------
// TODO: Perhaps I should rename this FESurfaceData (there is already a class called that though)
//       and then define FESurfaceMap as a tool for evaluating data across a surface (i.e. via shape functions)
class FESurfaceMap
{
public:
	//! default constructor
	FESurfaceMap(int dataType);

	//! copy constructor
	FESurfaceMap(const FESurfaceMap& map);

	//! assignment operator
	FESurfaceMap& operator = (const FESurfaceMap& map);

	//! Create a surface data map for this surface
	bool Create(const FESurface* ps);

	//! Create a surface data map for this surface
	bool Create(const FEFacetSet* ps);

	//! get the value for a given facet index
	template <class T> T GetValue(const FEFacetIndex& n) const;

	//! set the value
	template <class T> bool SetValue(const T& v);

	//! set the value
	template <class T> bool SetValue(const FEFacetIndex& n, const T& v);

	//! serialization
	void Serialize(DumpStream& ar);

	//! set the name
	void SetName(const std::string& name);

	//! get the name
	const std::string& GetName() const { return m_name; }

	//! get data type
	int DataType() const { return m_dataType; }

protected:
	int DataSize() const;

private:
	int	m_dataType;				//!< data type

	// default values
	double	m_defDouble;
	vec3d	m_defVec3d;

	std::vector<double>	m_val;	//!< data values
	std::string	m_name;
};

template <> inline double FESurfaceMap::GetValue<double>(const FEFacetIndex& n) const
{
	assert(m_dataType == FE_DOUBLE);
	return	m_val[n];
}

template <> inline vec3d FESurfaceMap::GetValue<vec3d>(const FEFacetIndex& n) const
{
	assert(m_dataType == FE_VEC3D);
	return	vec3d(m_val[3*n], m_val[3*n + 1], m_val[3*n+2]);
}
