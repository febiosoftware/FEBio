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
class FECORE_API FEDataGenerator : public FECoreBase
{
	FECORE_SUPER_CLASS

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

public:
	// overload  one of these functions for custom generators
	virtual void value(const vec3d& r, double& data) {}
	virtual void value(const vec3d& r, vec2d& data) {}
	virtual void value(const vec3d& r, vec3d& data) {}
	virtual void value(const vec3d& r, mat3d& data) {}
};
