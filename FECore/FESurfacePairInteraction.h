#pragma once
#include "FEModelComponent.h"
#include "FESurface.h"

//-----------------------------------------------------------------------------
class FEModel;
class FEGlobalMatrix;

//-----------------------------------------------------------------------------
//! This class describes a general purpose interaction between two surfaces.
class FESurfacePairInteraction : public FEModelComponent
{
public:
	//! constructor
	FESurfacePairInteraction(FEModel* pfem);

public:
	//! initialization routine
	virtual bool Init() = 0;

	//! update 
	virtual void Update(int niter) = 0;

	//! return the master surface
	virtual FESurface* GetMasterSurface() = 0;

	//! return the slave surface
	virtual FESurface* GetSlaveSurface () = 0;

	//! temporary construct to determine if contact interface uses nodal integration rule (or facet)
	virtual bool UseNodalIntegration() = 0;

	//! create a copy of this interface
	virtual void CopyFrom(FESurfacePairInteraction* pci) {}

	//! build the matrix profile for use in the stiffness matrix
	virtual void BuildMatrixProfile(FEGlobalMatrix& K) = 0;
};
