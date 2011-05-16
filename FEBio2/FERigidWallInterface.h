// FERigidWallInterface.h: interface for the FERigidWallInterface class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FERIGIDWALLINTERFACE_H__11D94285_BF60_4052_808A_562C563CA820__INCLUDED_)
#define AFX_FERIGIDWALLINTERFACE_H__11D94285_BF60_4052_808A_562C563CA820__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "FECore/FEContactInterface.h"
#include "FERigidSurface.h"
#include "FECore/FESurface.h"
#include "FECore/vec3d.h"
#include "FECore/vec2d.h"
#include "FECore/vector.h"

//-----------------------------------------------------------------------------
class FERigidWallSurface : public FESurface
{
public:
	//! constructor
	FERigidWallSurface(FEMesh* pm=0) : FESurface(pm) { m_NQ.Attach(this); }

	//! Initializes data structures
	void Init();

	//! shallow copy
	void ShallowCopy(FERigidWallSurface& s)
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

	void Serialize(DumpFile& ar);

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
//! This class implements a sliding contact interface with a rigid wall

//! This class is a specialization of the general sliding interface where
//! the master surface is a rigid wall

class FERigidWallInterface : public FEContactInterface
{
public:
	//! constructor
	FERigidWallInterface(FEModel* pfem);

	//! destructor
	virtual ~FERigidWallInterface() { if (m_mp) delete m_mp; }

	//! intializes rigid wall interface
	void Init();

	//! project slave nodes onto master plane
	void ProjectSurface(FERigidWallSurface& s);

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

	//! serialize data to archive
	void Serialize(DumpFile& ar);

private:
	//! copy constructor hidden
	FERigidWallInterface(FERigidWallInterface& ri){}

public:
	FERigidWallSurface	m_ss;		//!< slave surface
	FERigidSurface		*m_mp;		//!< master surface

	int nse;	//!< number of slave elements

	double	m_atol;	//!< augmentation tolerance

	double			m_eps;	//!< penalty scale factor
	int				m_nplc;	//!< penalty load curve number
	FELoadCurve*	m_pplc;	//!< pointer to penalty load curve

	DECLARE_PARAMETER_LIST();
};


#endif // !defined(AFX_FERIGIDWALLINTERFACE_H__11D94285_BF60_4052_808A_562C563CA820__INCLUDED_)
