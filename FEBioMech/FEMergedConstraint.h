#pragma once
#include <FECore/FEMesh.h>
#include "febiomech_api.h"
#include <vector>

//-----------------------------------------------------------------------------
// Forward declaration of FEModel class
class FEModel;

//-----------------------------------------------------------------------------
// A merged constraint takes two surfaces and merges them by matching each node of one surface
// to the corresponding node on the other surface and then generates a linear constraint 
// between the two nodes that essentially matches the degrees of freedom.
class FEBIOMECH_API FEMergedConstraint
{
public:
	FEMergedConstraint(FEModel& fem);
	~FEMergedConstraint();

	bool Merge(FEFacetSet* surf1, FEFacetSet* surf2, const std::vector<int>& dofList);

private:
	FEModel&	m_fem;
};
