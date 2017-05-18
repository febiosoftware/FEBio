// FESlidingInterface.h: interface for the FESlidingInterface class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FESLIDINGINTERFACE_H__742CFDED_B4BF_47AF_B673_1D26FE03F934__INCLUDED_)
#define AFX_FESLIDINGINTERFACE_H__742CFDED_B4BF_47AF_B673_1D26FE03F934__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "FEContactSurface.h"
#include "FEContactInterface.h"
#include "FECore/FEClosestPointProjection.h"
#include "FECore/vector.h"

//-----------------------------------------------------------------------------
class FESlidingSurface : public FEContactSurface
{
public:
	//! constructor
	FESlidingSurface(FEModel* pfem) : FEContactSurface(pfem) {}

	//! Initializes data structures
	bool Init();

	//! Calculate the total traction at a node
	vec3d traction(int inode);

	//! evaluate net contact force
	vec3d GetContactForce();
	
	//! evaluate net contact area
	double GetContactArea();
    
	//! Serialize data to archive
	void Serialize(DumpStream& ar);

public:
    void GetContactGap     (int nface, double& pg);
    void GetContactPressure(int nface, double& pg);
    void GetContactTraction(int nface, vec3d& pt);
	void GetNodalContactGap     (int nface, double* pg);
	void GetNodalContactPressure(int nface, double* pg);
	void GetNodalContactTraction(int nface, vec3d* pt);

public:
	vector<double>				m_gap;	//!< gap function at nodes
	vector<vec3d>				m_nu;	//!< master normal at slave node
	vector<FESurfaceElement*>	m_pme;	//!< master element a slave node penetrates
	vector<vec2d>				m_rs;	//!< natural coordinates of slave projection on master element
	vector<vec2d>				m_rsp;	//!< natural coordinates at previous time step
	vector<double>				m_Lm;	//!< Lagrange multipliers for contact pressure
	vector<mat2d>				m_M;	//!< surface metric tensor
	vector<vec2d>				m_Lt;	//!< Lagrange multipliers for friction
	vector<double>				m_off;	//!< gap offset (= shell thickness)
	vector<double>				m_eps;	//!< normal penalty factors
	vector<double>				m_Ln;	//!< net contact pressure
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
	FESlidingInterface(FEModel* pfem);

	//! destructor
	virtual ~FESlidingInterface(){}

	//! Initializes sliding interface
	bool Init();

	//! interface activation
	void Activate();

	//! update interface data
	void Update(int niter);

	//! projects slave nodes onto master nodes
	void ProjectSurface(FESlidingSurface& ss, FESlidingSurface& ms, bool bupseg, bool bmove = false);

	//! calculate penalty value
	double Penalty() { return m_eps; }

	//! calculate contact forces
	virtual void ContactForces(FEGlobalVector& R);

	//! calculate contact stiffness
	virtual void ContactStiffness(FESolver* psolver);

	//! calculate Lagrangian augmentations
	virtual bool Augment(int naug);

	//! serialize data to archive
	void Serialize(DumpStream& ar);

	//! calculate contact pressures for file output
	void UpdateContactPressures();

	//! return the master and slave surface
	FESurface* GetMasterSurface() { return &m_ms; }
	FESurface* GetSlaveSurface () { return &m_ss; }

	//! return integration rule class
	bool UseNodalIntegration() { return true; }

	//! build the matrix profile for use in the stiffness matrix
	void BuildMatrixProfile(FEGlobalMatrix& K);

protected:
	//! calculate auto penalty factor
	void CalcAutoPenalty(FESlidingSurface& s);

	//! calculate the nodal force of a slave node
	void ContactNodalForce(int m, FESlidingSurface& ss, FESurfaceElement& mel, vector<double>& fe);

	//! calculate the stiffness contribution of a single slave node
	void ContactNodalStiffness(int m, FESlidingSurface& ss, FESurfaceElement& mel, matrix& ke);

	//! map the frictional data from the old element to the new element
	void MapFrictionData(int inode, FESlidingSurface& ss, FESlidingSurface& ms, FESurfaceElement& sn, FESurfaceElement& so, vec3d& q);

public:
	FESlidingSurface	m_ss;	//!< slave surface
	FESlidingSurface	m_ms;	//!< master surface

	bool			m_btwo_pass;	//!< two pass algorithm flag

	bool			m_bself_contact;	//!< self-contact flag
	double			m_sradius;			//!< search radius for self contact

	int				m_naugmax;	//!< maximum nr of augmentations
	int				m_naugmin;	//!< minimum nr of augmentations
	double			m_gtol;		//!< gap tolerance
	double			m_atol;		//!< augmentation tolerance
	double			m_ktmult;	//!< multiplier for tangential stiffness
	double			m_knmult;	//!< multiplier for normal stiffness

	double			m_stol;		//!< search tolerance

	bool			m_bautopen;	//!< auto penalty calculation
	double			m_eps;		//!< penalty scale factor 

	bool			m_breloc;	//!< initial node relocation

	double			m_mu;		//!< friction coefficient
	double			m_epsf;		//!< penalty scale factor for friction

	int				m_nsegup;	//!< segment update parameter

private:
	bool	m_bfirst;	//!< flag to indicate the first time we enter Update
	double	m_normg0;	//!< initial gap norm

public:
	DECLARE_PARAMETER_LIST();
};

#endif // !defined(AFX_FESLIDINGINTERFACE_H__742CFDED_B4BF_47AF_B673_1D26FE03F934__INCLUDED_)
