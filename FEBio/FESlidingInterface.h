// FESlidingInterface.h: interface for the FESlidingInterface class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FESLIDINGINTERFACE_H__742CFDED_B4BF_47AF_B673_1D26FE03F934__INCLUDED_)
#define AFX_FESLIDINGINTERFACE_H__742CFDED_B4BF_47AF_B673_1D26FE03F934__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "FEContactInterface.h"
#include "FESurface.h"
#include "vec2d.h"

//-----------------------------------------------------------------------------
class FESlidingSurface : public FESurface
{
public:
	//! constructor
	FESlidingSurface(FEMesh* pm=0) : FESurface(pm) { m_NQ.Attach(this); }

	//! Initializes data structures
	void Init();

	//! shallow copy
	void ShallowCopy(FESlidingSurface& s)
	{
		Lm  = s.Lm;
		gap = s.gap;
		zero(pme);
		Lt  = s.Lt;
	}

	//! Update the surface data
	void Update() {}

	//! Find element that contains the projection of x
	FEElement* FindMasterSegment(vec3d& x, vec3d& q, vec2d& r, bool& binit_nq, double tol);

	//! Calculate the total traction at a node
	vec3d traction(int inode);

	void UpdateNormals();

public:
	vector<double>		gap;	//!< gap function at nodes
	vector<vec3d>		nu;		//!< master normal at slave node
	vector<FEElement*>	pme;	//!< master element a slave node penetrates
	vector<vec2d>		rs;		//!< natural coordinates of slave projection on master element
	vector<vec2d>		rsp;	//!< natural coordinates at previous time step
	vector<double>		Lm;		//!< Lagrange multipliers for contact pressure
	vector<mat2d>		M;		//!< surface metric tensor
	vector<vec2d>		Lt;		//!< Lagrange multipliers for friction
	vector<double>		off;	//!< gap offset (= shell thickness)
	vector<double>		eps;	//!< normal penalty factors

	FENNQuery		m_NQ;		//!< this structure is used in finding the master element that corresponds to a slave node
};

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
	void ProjectSurface(FESlidingSurface& ss, FESlidingSurface& ms, bool bupseg, bool bmove = false);

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

	//! serialize data to archive
	void Serialize(Archive& ar);

protected:
	//! calculate auto penalty factor
	void CalcAutoPenalty(FESlidingSurface& s);

	//! calculate the nodal force of a slave node
	void ContactNodalForce(int m, FESlidingSurface& ss, FESurfaceElement& mel, vector<double>& fe);

	//! calculate the stiffness contribution of a single slave node
	void ContactNodalStiffness(int m, FESlidingSurface& ss, FESurfaceElement& mel, matrix& ke);

	//! map the frictional data from the old element to the new element
	void MapFrictionData(int inode, FESlidingSurface& ss, FESlidingSurface& ms, FESurfaceElement& sn, FESurfaceElement& so, vec3d& q);

private:
	//! copy constructor hidden
	FESlidingInterface(FESlidingInterface& si){}

public:
	FESlidingSurface	m_ss;	//!< slave surface
	FESlidingSurface	m_ms;	//!< master surface

	int				m_npass;	//!< nr of passes (1 or 2)

	int				m_naugmax;	//!< maximum nr of augmentations
	int				m_naugmin;	//!< minimum nr of augmentations
	double			m_gtol;		//!< gap tolerance
	double			m_atol;		//!< augmentation tolerance
	double			m_ktmult;	//!< multiplier for tangential stiffness
	double			m_knmult;	//!< multiplier for normal stiffness

	double			m_stol;		//!< search tolerance

	int				m_nautopen;	//!< auto penalty calculation factor (0=none, 1=old, 2=new)
	double			m_eps;		//!< penalty scale factor 

	bool			m_breloc;	//!< initial node relocation

	double			m_mu;		//!< friction coefficient
	double			m_epsf;		//!< penalty scale factor for friction

	int				m_nsegup;	//!< segment update parameter
	int				m_nplc;		//!< penalty load curve number
	FELoadCurve*	m_pplc;		//!< pointer to penalty load curve
};

#endif // !defined(AFX_FESLIDINGINTERFACE_H__742CFDED_B4BF_47AF_B673_1D26FE03F934__INCLUDED_)
