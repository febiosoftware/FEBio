#pragma once
#include "FEBoundaryCondition.h"

//-----------------------------------------------------------------------------
class FENodeSet;
class FEFacetSet;

//-----------------------------------------------------------------------------
// base class for prescribed boundary conditions
class FECORE_API FEPrescribedBC : public FEBoundaryCondition
{
public:
	FEPrescribedBC(FEModel* pfem);

public:
	// implement these functions

	// assign a node set to the prescribed BC
	virtual void AddNodes(const FENodeSet& set) {};

	// assign a surface to the BC
	// By default, the nodes of the surface are assigned to the BC
	virtual void AddNodes(const FEFacetSet& surf);

	// This function is called when the solver needs to know the 
	// prescribed dof values. The brel flag indicates wheter the total 
	// value is needed or the value with respect to the current nodal dof value
	virtual void PrepStep(std::vector<double>& ui, bool brel = true) = 0;

	// copy data from another class
	virtual void CopyFrom(FEPrescribedBC* pbc) = 0;
};
