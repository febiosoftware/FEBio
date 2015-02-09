#pragma once
#include "FEModelComponent.h"
#include "FESurface.h"

//-----------------------------------------------------------------------------
class FEModel;

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

	//! Create a shallow copy
	virtual void ShallowCopy(DumpStream& dmp, bool bsave) = 0;

	//! serialize data to archive
	virtual void Serialize(DumpFile& ar) = 0;

	//! return the master surface
	virtual FESurface* GetMasterSurface() = 0;

	//! return the slave surface
	virtual FESurface* GetSlaveSurface () = 0;

	//! temporary construct to determine if contact interface uses nodal integration rule (or facet)
	virtual bool UseNodalIntegration() = 0;

	//! create a copy of this interface
	virtual void CopyFrom(FESurfacePairInteraction* pci) {}

	//! Get the Contact Interface ID
	int GetID() { return m_nID; }

protected:
	int		m_nID;			//!< ID of interface

	static int	m_ncount;	//!< used to create unique ID's for the contact interfaces
};
