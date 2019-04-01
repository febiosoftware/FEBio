#pragma once
#include "FEElasticMaterial.h"
#include "FEElasticFiberMaterial.h"
#include "FEFiberDensityDistribution.h"

//----------------------------------------------------------------------------------
// This is an iterator class that can be used to loop over all integration points of
// a fiber integration scheme. It must implement the Next functions.
// During the Next call, it should also update the fiber vector and weight.
// Next should return false if there is no more integration point
class FEFiberIntegrationSchemeIterator
{
public:
	FEFiberIntegrationSchemeIterator() {}
	virtual ~FEFiberIntegrationSchemeIterator() {}

	// Move to the next integration point
	// This also updates the m_fiber and m_weight members
	virtual bool Next() = 0;

	// check if the iterator is valid
	virtual bool IsValid() = 0;

public:
	vec3d	m_fiber;		// current fiber vector at integration point
	double	m_weight;		// current integration weight
};

//----------------------------------------------------------------------------------
// Base clase for integration schemes for continuous fiber distributions.
// The purpose of this class is mainly to provide an interface to the integration schemes
// for the FEBio input file. The code will use the GetIterator function to create an
// iterator that can be used to loop over all the integration points of the scheme and to
// evaluate the fiber vector and weights at each point.
class FEFiberIntegrationScheme : public FEMaterial
{
public:
    FEFiberIntegrationScheme(FEModel* pfem);
    
	// Creates an iterator for the scheme. 
	// In general, the integration scheme may depend on the material point.
	// The passed material point pointer will be zero when evaluating the integrated fiber density
	virtual FEFiberIntegrationSchemeIterator* GetIterator(FEMaterialPoint* mp = 0) = 0;
};
