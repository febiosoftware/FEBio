// FERigidWallInterface.h: interface for the FERigidWallInterface class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FERIGIDWALLINTERFACE_H__11D94285_BF60_4052_808A_562C563CA820__INCLUDED_)
#define AFX_FERIGIDWALLINTERFACE_H__11D94285_BF60_4052_808A_562C563CA820__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "FEContactInterface.h"
#include "FERigidSurface.h"
#include "FEContactSurface.h"

//-----------------------------------------------------------------------------
//! This class implements a sliding contact interface with a rigid wall

//! This class is a specialization of the general sliding interface where
//! the master surface is a rigid wall

class FERigidWallInterface : public FEContactInterface
{
public:
	//! constructor
	FERigidWallInterface(FEM* pfem);

	//! destructor
	virtual ~FERigidWallInterface() { if (m_mp) delete m_mp; }

	//! intializes rigid wall interface
	void Init();

	//! project slave nodes onto master plane
	void ProjectSurface(FEContactSurface& s);

	//! update rigid wall data
	void Update();

	//! shallow copy
	void ShallowCopy(FEContactInterface& ci)
	{
		FERigidWallInterface& ri = dynamic_cast<FERigidWallInterface&>(ci);
		m_ss.ShallowCopy(ri.m_ss);
	}

	//! calculate penalty value
	double Penalty() { return (m_pplc?m_eps*m_pplc->Value():m_eps);}

	//! calculate contact forces
	virtual void ContactForces(vector<double>& F);

	//! calculate contact stiffness
	virtual void ContactStiffness();

	//! calculate Lagrangian augmentations
	virtual bool Augment(int naug);

	//! Set the master surface
	void SetMasterSurface(FERigidSurface* prs) { assert(m_mp == 0); m_mp = prs; }

private:
	//! copy constructor hidden
	FERigidWallInterface(FERigidWallInterface& ri){}

public:
	FEContactSurface	m_ss;		//!< slave surface
	FERigidSurface		*m_mp;		//!< master surface

	int nse;	//!< number of slave elements

	double	m_atol;	//!< augmentation tolerance

	double			m_eps;	//!< penalty scale factor
	int				m_nplc;	//!< penalty load curve number
	FELoadCurve*	m_pplc;	//!< pointer to penalty load curve
};


#endif // !defined(AFX_FERIGIDWALLINTERFACE_H__11D94285_BF60_4052_808A_562C563CA820__INCLUDED_)
