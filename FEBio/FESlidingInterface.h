// FESlidingInterface.h: interface for the FESlidingInterface class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FESLIDINGINTERFACE_H__742CFDED_B4BF_47AF_B673_1D26FE03F934__INCLUDED_)
#define AFX_FESLIDINGINTERFACE_H__742CFDED_B4BF_47AF_B673_1D26FE03F934__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "FEContactInterface.h"
#include "FEContactSurface.h"

//-----------------------------------------------------------------------------
//! This class implements a sliding interface

//! The FESlidingInterface class defines an interface between two surfaces.
//! These surfaces define the slave and master surfaces of the contact
//! interface. 

class FESlidingInterface : public FEContactInterface 
{
public:
	//! constructor
	FESlidingInterface(FEM* pfem);

	//! destructor
	virtual ~FESlidingInterface(){}

	//! Initializes sliding interface
	void Init();

	//! update interface data
	void Update();

	//! projects slave nodes onto master nodes
	void ProjectSurface(FEContactSurface& ss, FEContactSurface& ms, bool bupseg, bool bmove = false);

	//! shallow copy
	void ShallowCopy(FEContactInterface& ci)
	{
		FESlidingInterface& si = dynamic_cast<FESlidingInterface&>(ci);
		m_ss.ShallowCopy(si.m_ss);
		m_ms.ShallowCopy(si.m_ms);
	}

	//! calculate penalty value
	double Penalty() { return (m_pplc?m_eps*m_pplc->Value():m_eps);}

	//! calculate contact forces
	virtual void ContactForces(vector<double>& F);

	//! calculate contact stiffness
	virtual void ContactStiffness();

	//! calculate Lagrangian augmentations
	virtual bool Augment(int naug);

protected:
	//! calculate auto penalty factor
	void CalcAutoPenalty(FEContactSurface& s);

	//! calculate the nodal force of a slave node
	void ContactNodalForce(int m, FEContactSurface& ss, FESurfaceElement& mel, vector<double>& fe);

	//! calculate the stiffness contribution of a single slave node
	void ContactNodalStiffness(int m, FEContactSurface& ss, FESurfaceElement& mel, matrix& ke);

	//! map the frictional data from the old element to the new element
//	void MapFrictionData(int inode, FESurfaceElement& sn, FESurfaceElement& so, vec3d& q);

private:
	//! copy constructor hidden
	FESlidingInterface(FESlidingInterface& si){}

public:
	FEContactSurface	m_ss;	//!< slave surface
	FEContactSurface	m_ms;	//!< master surface

	int				m_npass;	//!< nr of passes (1 or 2)

	int				m_naugmax;	//!< maximum nr of augmentations
	int				m_naugmin;	//!< minimum nr of augmentations
	double			m_gtol;		//!< gap tolerance
	double			m_atol;		//!< augmentation tolerance
	bool			m_blaugon;	//!< augmented lagrangian flag

	double			m_stol;		//!< search tolerance

	bool			m_bautopen;	//!< auto penalty calculation factor
	double			m_eps;		//!< penalty scale factor

	double			m_mu;		//!< friction coefficient
	double			m_epsf;		//!< penalty scale factor for friction

	int				m_nsegup;	//!< segment update parameter
	int				m_nplc;		//!< penalty load curve number
	FELoadCurve*	m_pplc;		//!< pointer to penalty load curve
};



#endif // !defined(AFX_FESLIDINGINTERFACE_H__742CFDED_B4BF_47AF_B673_1D26FE03F934__INCLUDED_)
