// FETiedInterface.h: interface for the FETiedInterface class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FETIEDINTERFACE_H__2D0F2799_9B2D_463F_B42E_C6924D6BCD6E__INCLUDED_)
#define AFX_FETIEDINTERFACE_H__2D0F2799_9B2D_463F_B42E_C6924D6BCD6E__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "FEContactInterface.h"
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
	bool Init() override;

	//! interface activation
	void Activate() override;

	//! projects slave nodes onto master nodes
	void ProjectSurface(FETiedContactSurface& ss, FETiedContactSurface& ms, bool bmove = false);

	//! serialize data to archive
	void Serialize(DumpStream& ar) override;

	//! return the master and slave surface
	FESurface* GetMasterSurface() override { return &ms; }
	FESurface* GetSlaveSurface () override { return &ss; }

	//! return integration rule class
	bool UseNodalIntegration() override { return true; }

	//! build the matrix profile for use in the stiffness matrix
	void BuildMatrixProfile(FEGlobalMatrix& K) override;

public:
	//! calculate contact forces
	void Residual(FEGlobalVector& R, const FETimeInfo& tp) override;

	//! calculate contact stiffness
	void StiffnessMatrix(FESolver* psolver, const FETimeInfo& tp) override;

	//! calculate Lagrangian augmentations
	bool Augment(int naug, const FETimeInfo& tp) override;

	//! update
	void Update() override;

public:
	FETiedContactSurface	ss;	//!< slave surface
	FETiedContactSurface	ms;	//!< master surface

public:
	double		m_atol;		//!< augmentation tolerance
	double		m_eps;		//!< penalty scale factor
	double		m_stol;		//!< search tolerance
	int			m_naugmax;	//!< maximum nr of augmentations
	int			m_naugmin;	//!< minimum nr of augmentations
	bool		m_boffset;	//!< offset slave surface for shells
	double		m_Dmax;		//!< max distance for contact
	bool		m_bspecial;	//!< handle special cases in projection
	bool		m_breloc;	//!< node relocation on initialization

	DECLARE_FECORE_CLASS();
};

#endif // !defined(AFX_FETIEDINTERFACE_H__2D0F2799_9B2D_463F_B42E_C6924D6BCD6E__INCLUDED_)
