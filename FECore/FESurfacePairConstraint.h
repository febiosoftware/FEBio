#pragma once
#include "FEModelComponent.h"
#include "FESurface.h"

//-----------------------------------------------------------------------------
class FEModel;
class FEGlobalMatrix;

//-----------------------------------------------------------------------------
//! This class describes a general purpose interaction between two surfaces.
// TODO: I like to inherit this from FENLConstraint and potentially eliminate
//       The distinction between a nonlinear constraint and a contact interface.
//       Since a contact interface essentially is a nonlinear constraint, I think
//       this may make things a lot easier. I already made the function definitions consistent
//       but am hesitant to push this through at this point. 
class FECORE_API FESurfacePairConstraint : public FEModelComponent
{
	DECLARE_SUPER_CLASS(FESURFACEPAIRINTERACTION_ID);

public:
	//! constructor
	FESurfacePairConstraint(FEModel* pfem);

public:
	//! return the master surface
	virtual FESurface* GetMasterSurface() = 0;

	//! return the slave surface
	virtual FESurface* GetSlaveSurface () = 0;

	//! temporary construct to determine if contact interface uses nodal integration rule (or facet)
	virtual bool UseNodalIntegration() = 0;

	//! create a copy of this interface
	virtual void CopyFrom(FESurfacePairConstraint* pci) {}

public:
	// Build the matrix profile
	virtual void BuildMatrixProfile(FEGlobalMatrix& M) = 0;

	// reset the state data
	virtual void Reset() {}
};
