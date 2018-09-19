#pragma once
#include "FENodeDataMap.h"
#include "FESurfaceMap.h"
#include "FEDomainMap.h"
#include "FEElementSet.h"
#include <string>

class FENodeSet;
class FEFacetSet;

//-----------------------------------------------------------------------------
class FECORE_API FEDataMathGenerator
{
public:
	FEDataMathGenerator();

	// set the math expression
	void setExpression(const std::string& math);

	// generate the data array for the given node set
	bool Generate(FENodeDataMap& ar, const FENodeSet& set);

	// generate the data array for the given facet set
	bool Generate(FESurfaceMap& data, const FEFacetSet& surf);

	// generate the data array for the given element set
	bool Generate(FEDomainMap& data, FEElementSet& set);

private:
	std::string	m_math;
};
