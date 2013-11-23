#pragma once
#include "FEParameterList.h"
#include "FESurface.h"

//-----------------------------------------------------------------------------
class FEModel;

//-----------------------------------------------------------------------------
//! This class describes a general purpose interaction between two surfaces.
class FESurfacePairInteraction : public FEParamContainer
{
public:
	//! constructor
	FESurfacePairInteraction(FEModel* pfem);

	//! return the type of this interface
	int Type();

	//! Get the FE model
	FEModel* GetFEModel() { return m_pfem; }

public:
	//! initialization routine
	virtual bool Init() = 0;

	//! update 
	virtual void Update(int niter) = 0;

	//! Create a shallow copy
	virtual void ShallowCopy(FESurfacePairInteraction& ci) = 0;

	//! serialize data to archive
	virtual void Serialize(DumpFile& ar) = 0;

	//! return the master surface
	virtual FESurface* GetMasterSurface() = 0;

	//! return the slave surface
	virtual FESurface* GetSlaveSurface () = 0;

public:
	//! Is this contact interface active
	bool IsActive();

	//! Activate the contact interface
	virtual void Activate();

	//! Deactivate the contact interface
	virtual void Deactivate();

protected:
	FEModel*	m_pfem;		//!< FEModel class this interface belongs to

protected:
	int		m_ntype;		//!< type of interface
	int		m_nID;			//!< ID of interface
	bool	m_bactive;		//!< active flag
};
