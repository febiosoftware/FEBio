#pragma once
#include "vec3d.h"
#include "FECoreBase.h"

//-----------------------------------------------------------------------------
class FENodeSet;
class FEFacetSet;
class FEElementSet;
class FENodeDataMap;
class FESurfaceMap;
class FEDomainMap;

//-----------------------------------------------------------------------------
// Data generators are used to generate values of model parameters. 
class FEDataGenerator : public FECoreBase
{
public:
	FEDataGenerator(FEModel* fem);
	virtual ~FEDataGenerator();

	// this function gives the data generator a chance to initialize itself
	// and check for any input problems.
	virtual bool Init();

	// generate the data array for the given node set
	bool Generate(FENodeDataMap& ar, const FENodeSet& set);

	// generate the data array for the given facet set
	bool Generate(FESurfaceMap& data, const FEFacetSet& surf);

	// generate the data array for the given element set
	bool Generate(FEDomainMap& data, FEElementSet& set);

	// get the FE model
	FEModel* GetFEModel() { return m_fem; }

public:
	// overload  this function for custom generators
	virtual void value(const vec3d& r, vector<double>& data) = 0;

private:
	FEModel*	m_fem;
};
