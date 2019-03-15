#pragma once
#include "FECoreBase.h"

//-----------------------------------------------------------------------------
// forward declarations
class FEModel;
class FEElement;

//-----------------------------------------------------------------------------
// Base class for all mesh adaptors
class FECORE_API FEMeshAdaptor : public FECoreBase
{
	FECORE_SUPER_CLASS

public:
	FEMeshAdaptor(FEModel* fem);

	// The mesh adaptor should return true if the mesh remained unchanged
	// otherwise, it should return false.
	// iteration is the iteration number of the mesh adaptation loop
	virtual bool Apply(int iteration) = 0;
};

//-----------------------------------------------------------------------------
// This class is a helper class for use in the mesh adaptors. Its purpose is to select
// elements based on some criterion. This element list is then usually passed to the 
// mesh adaptor.
class FECORE_API FEMeshAdaptorCriterion : public FECoreBase
{
	FECORE_SUPER_CLASS

public:
	FEMeshAdaptorCriterion(FEModel* fem);

	void SetSort(bool b);

	void SetMaxElements(int m);

public:

	// return a list of elements that satisfy the criterion
	virtual std::vector<int> GetElementList();

	// This function needs to be overridden in order to select some elements
	// that satisfy the selection criterion
	// return true if the element satisfies the criterion, otherwise false
	// If this function returns true, the elemVal paramter should be set
	// This is used to sort the element list
	virtual bool Check(FEElement& el, double& elemVal);

private:
	bool	m_sortList;		// sort the list
	int		m_maxelem;		// the max nr of elements to return (or 0 if don't care)

	DECLARE_FECORE_CLASS();
};
