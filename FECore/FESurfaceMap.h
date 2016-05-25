#pragma once
#include <vector>
#include <string>
#include <FECore/FEDataArray.h>
#include <assert.h>

//-----------------------------------------------------------------------------
class FESurface;
class FEFacetSet;
class DumpStream;

//-----------------------------------------------------------------------------
typedef int FEFacetIndex;


//-----------------------------------------------------------------------------
// TODO: Perhaps I should rename this FESurfaceData (there is already a class called that though)
//       and then define FESurfaceMap as a tool for evaluating data across a surface (i.e. via shape functions)
class FESurfaceMap : public FEDataArray
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

	//! serialization
	void Serialize(DumpStream& ar);

	//! set the name
	void SetName(const std::string& name);

	//! get the name
	const std::string& GetName() const { return m_name; }

private:
	std::string	m_name;
};
