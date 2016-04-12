// FERigidWallInterface.h: interface for the FERigidWallInterface class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FERIGIDWALLINTERFACE_H__11D94285_BF60_4052_808A_562C563CA820__INCLUDED_)
#define AFX_FERIGIDWALLINTERFACE_H__11D94285_BF60_4052_808A_562C563CA820__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "FECore/vector.h"
#include "FEContactInterface.h"
#include "FECore/FESurface.h"
#include "FECore/vec3d.h"
#include "FECore/vec2d.h"
#include "FECore/FENNQuery.h"
#include "FERigidSurface.h"

//-----------------------------------------------------------------------------
class FERigidWallSurface : public FESurface
{
public:
	//! constructor
	FERigidWallSurface(FEModel* pfem);

	//! Initializes data structures
	bool Init();

	//! Update the surface data
	void Update() {}

	//! Calculate the total traction at a node
	vec3d traction(int inode);

	void UpdateNormals();

	void Serialize(DumpStream& ar);

	//! Unpack surface element data
	void UnpackLM(FEElement& el, vector<int>& lm);

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

	int	m_dofX;
	int	m_dofY;
	int	m_dofZ;

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
	bool Init();

	//! interface activation
	void Activate();

	//! project slave nodes onto master plane
	void ProjectSurface(FERigidWallSurface& s);

	//! update rigid wall data
	void Update(int niter);

	//! calculate contact forces
	virtual void ContactForces(FEGlobalVector& R);

	//! calculate contact stiffness
	virtual void ContactStiffness(FESolver* psolver);

	//! calculate Lagrangian augmentations
	virtual bool Augment(int naug);

	//! Set the master surface
	void SetMasterSurface(FERigidSurface* prs) { assert(m_mp == 0); m_mp = prs; }

	//! serialize data to archive
	void Serialize(DumpStream& ar);

	//! return the master and slave surface
	FESurface* GetMasterSurface() { return 0; }
	FESurface* GetSlaveSurface () { return &m_ss; }

	//! return integration rule class
	bool UseNodalIntegration() { return true; }

	//! build the matrix profile for use in the stiffness matrix
	void BuildMatrixProfile(FEGlobalMatrix& K);

public: // inherited from FECoreBase

	//! get the number of properties
	int Properties();

	//! get a specific property
	FECoreBase* GetProperty(int i);

	//! find a property index ( returns <0 for error)
	int FindPropertyIndex(const char* szname);

	//! set a property (returns false on error)
	bool SetProperty(int i, FECoreBase* pm);

public:
	FERigidWallSurface	m_ss;		//!< slave surface
	FERigidSurface		*m_mp;		//!< master surface

	int nse;	//!< number of slave elements

	double		m_atol;		//!< augmentation tolerance
	double		m_eps;		//!< penalty scale factor
	double		m_d;		//!< normal offset

	DECLARE_PARAMETER_LIST();
};


#endif // !defined(AFX_FERIGIDWALLINTERFACE_H__11D94285_BF60_4052_808A_562C563CA820__INCLUDED_)
