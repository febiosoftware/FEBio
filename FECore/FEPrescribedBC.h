#pragma once
#include "FEBoundaryCondition.h"

//-----------------------------------------------------------------------------
class FENodeSet;
class FESurface;

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
	virtual void AddNodes(const FESurface& surf) {}

	// This function is called when the solver needs to know the 
	// prescribed dof values. The brel flag indicates wheter the total 
	// value is needed or the value with respect to the current nodal dof value
	virtual void PrepStep(std::vector<double>& ui, bool brel = true) = 0;

	// This is called during nodal update and should be used to enforce the 
	// nodal degrees of freedoms
	virtual void Update() = 0;

	// copy data from another class
	virtual void CopyFrom(FEPrescribedBC* pbc) = 0;

	// Also implement the following functions.
	// These are already declared in base classes.
	//  bool Init();
	//  void Activate();
	//  void Deactivate();
	//  void Serialize(DumpStream& ar);
};
