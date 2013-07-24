// FETiedInterface.h: interface for the FETiedInterface class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FETIEDINTERFACE_H__2D0F2799_9B2D_463F_B42E_C6924D6BCD6E__INCLUDED_)
#define AFX_FETIEDINTERFACE_H__2D0F2799_9B2D_463F_B42E_C6924D6BCD6E__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "FECore/FEContactInterface.h"
#include "FETiedContactSurface.h"

//-----------------------------------------------------------------------------
//! This class implements a tied interface.

class FETiedInterface : public FEContactInterface
{
public:
	//! constructor
	FETiedInterface(FEModel* pfem);

	//! destructor
	virtual ~FETiedInterface(){}

	//! Initializes sliding interface
	bool Init();

	//! interface activation
	void Activate();

	//! update interface data
	void Update(int niter);

	//! projects slave nodes onto master nodes
	void ProjectSurface(FETiedContactSurface& ss, FETiedContactSurface& ms, bool bmove = false);

	//! shallow copy
	void ShallowCopy(FEContactInterface& ci);

	//! calculate contact forces
	virtual void ContactForces(FEGlobalVector& R);

	//! calculate contact stiffness
	virtual void ContactStiffness(FESolver* psolver);

	//! calculate Lagrangian augmentations
	virtual bool Augment(int naug);

	//! serialize data to archive
	void Serialize(DumpFile& ar);

	//! return the master and slave surface
	FESurface* GetMasterSurface() { return &ms; }
	FESurface* GetSlaveSurface () { return &ss; }

private:
	//! copy constructor hidden
	FETiedInterface(FETiedInterface& si){}

public:
	FETiedContactSurface	ss;	//!< slave surface
	FETiedContactSurface	ms;	//!< master surface

	int nse;	//!< number of slave elements
	int nme;	//!< number of master elements

	double		m_atol;		//!< augmentation tolerance
	double		m_eps;		//!< penalty scale factor
	double		m_stol;		//!< search tolerance
	int			m_naugmax;	//!< maximum nr of augmentations
	int			m_naugmin;	//!< minimum nr of augmentations

	DECLARE_PARAMETER_LIST();
};

#endif // !defined(AFX_FETIEDINTERFACE_H__2D0F2799_9B2D_463F_B42E_C6924D6BCD6E__INCLUDED_)
