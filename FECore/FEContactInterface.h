// FEContactInterface.h: interface for the FEContactInterface class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FECONTACTINTERFACE_H__DB89AB6C_87DA_4C47_923C_B9D488E04403__INCLUDED_)
#define AFX_FECONTACTINTERFACE_H__DB89AB6C_87DA_4C47_923C_B9D488E04403__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "FEMesh.h"
#include "DumpFile.h"
#include "FEParameterList.h"
#include "FESolver.h"

class FEModel;

// Macauley bracket
#define MBRACKET(x) ((x)>=0? (x): 0)

// Heavyside function
#define HEAVYSIDE(x) ((x)>=0?1:0)


//-----------------------------------------------------------------------------
//! This is the base class for contact interfaces
//! \todo Merge the active interface with the boundary condition class
class FEContactInterface : public FEParamContainer
{
public:
	//! constructor
	FEContactInterface(FEModel* pfem);

	//! destructor
	virtual ~FEContactInterface();

	//! initialization routine
	virtual bool Init() = 0;

	//! update 
	virtual void Update(int niter) = 0;

	//! return the type of this interface
	int Type() { return m_ntype; }

	//! Create a shallow copy
	virtual void ShallowCopy(FEContactInterface& ci)=0;

	//! calculate contact forces
	virtual void ContactForces(FEGlobalVector& R) = 0;

	//! calculate contact stiffness
	virtual void ContactStiffness(FESolver* psolver) = 0;

	//! calculate Lagrangian augmentations
	virtual bool Augment(int naug) = 0;

	//! serialize data to archive
	virtual void Serialize(DumpFile& ar);

	//! return the master and slave surface
	virtual FESurface* GetMasterSurface() = 0;
	virtual FESurface* GetSlaveSurface () = 0;

public:
	//! Is this contact interface active
	bool IsActive() { return m_bactive; }

	//! Activate the contact interface
	virtual void Activate() { m_bactive = true; }

	//! Deactivate the contact interface
	virtual void Deactivate() { m_bactive = false; }

protected:
	//! don't call the default constructor
	FEContactInterface() {m_pfem=0;}

	//! auto-penalty calculation
	double AutoPenalty(FESurfaceElement& el, FESurface& s);

public:
	bool	m_blaugon;	//!< augmented lagrangian flag

protected:
	FEModel*	m_pfem;	//!< FEModel class this interface belongs to

	int	m_ntype;	//!< type of interface

	int	m_nID;	//!< ID of interface

protected:
	bool	m_bactive;	//!< active flag
};

#endif // !defined(AFX_FECONTACTINTERFACE_H__DB89AB6C_87DA_4C47_923C_B9D488E04403__INCLUDED_)
