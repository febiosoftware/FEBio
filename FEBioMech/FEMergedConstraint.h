#pragma once
#include <FECore/FEMesh.h>
#include <vector>
using namespace std;

class FEModel;

//-----------------------------------------------------------------------------
// A merged constraint takes two surfaces and merges them by matching each node of one surface
// to the corresponding node on the other surface and then generates a linear constraint 
// between the two nodes that essentially matches the degrees of freedom.
class FECORE_API FEMergedConstraint
{
public:
	FEMergedConstraint(FEModel& fem);
	~FEMergedConstraint();

	bool Merge(FEFacetSet* surf1, FEFacetSet* surf2, const vector<int>& dofList);

private:
	FEModel&	m_fem;
};
