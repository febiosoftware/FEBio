#pragma once
#include <vector>

//-----------------------------------------------------------------------------
class FESurface;
class DumpStream;

//-----------------------------------------------------------------------------
typedef int FEFacetIndex;

//-----------------------------------------------------------------------------
// TODO: Perhaps I should rename this FESurfaceData (there is already a class called that though)
//       and then define FESurfaceMap as a tool for evaluating data across a surface (i.e. via shape functions)
class FESurfaceMap
{
public:
	FESurfaceMap();

	//! Create a surface data map for this surface
	bool Create(const FESurface* ps, double def = 0.0);

	//! get the value for a given facet index
	double GetValue(const FEFacetIndex& n) const { return m_val[n]; }

	//! set the value
	bool SetValue(const FEFacetIndex& n, double v);

	//! serialization
	void Serialize(DumpStream& ar);

private:
	std::vector<double>	m_val;	//!< data values
};
