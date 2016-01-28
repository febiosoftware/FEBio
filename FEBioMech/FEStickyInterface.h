#pragma once
#include "FEContactInterface.h"
#include "FEContactSurface.h"

//-----------------------------------------------------------------------------
//! This class describes a contact slave or master surface used for 
//! sticky contact

//!	this class is used in contact analyses to describe a contacting
//! surface in a sticky contact interface.

class FEStickySurface : public FEContactSurface
{
public:
	class NODE 
	{
	public:
		NODE() { gap = vec3d(0.0,0.0,0.0); pme = 0; }

	public:
		vec3d				gap;	//!< "gap" function
		vec2d				rs;		//!< natural coordinates of slave projection on master element
		vec3d				Lm;		//!< Lagrange multiplier
		vec3d				tn;		//!< traction vector
		FESurfaceElement*	pme;	//!< master element a slave node penetrates
	};

public:
	//! constructor
	FEStickySurface(FEModel* pfem) : FEContactSurface(pfem) {}

	//! Initializes data structures
	bool Init();

	//! data serialization
	void Serialize(DumpStream& ar);

public:
	void GetNodalContactGap     (int nface, double* gn);
	void GetNodalContactPressure(int nface, double* pn);
	void GetNodalContactTraction(int nface, vec3d* tn);

public:
	vector<NODE>	m_Node;	//!< node contact data
};

//-----------------------------------------------------------------------------
//! This class implements a sticky interface.
//! A sticky interface is like tied, but nodes are only tied when they come into
//! contact.

class FEStickyInterface : public FEContactInterface
{
public:
	//! constructor
	FEStickyInterface(FEModel* pfem);

	//! destructor
	virtual ~FEStickyInterface(){}

	//! Initializes sliding interface
	bool Init();

	//! interface activation
	void Activate();

	//! update interface data
	void Update(int niter);

	//! projects slave nodes onto master nodes
	void ProjectSurface(FEStickySurface& ss, FEStickySurface& ms, bool bmove = false);

	//! calculate contact forces
	virtual void ContactForces(FEGlobalVector& R);

	//! calculate contact stiffness
	virtual void ContactStiffness(FESolver* psolver);

	//! calculate Lagrangian augmentations
	virtual bool Augment(int naug);

	//! serialize data to archive
	void Serialize(DumpStream& ar);

	//! return the master and slave surface
	FESurface* GetMasterSurface() { return &ms; }
	FESurface* GetSlaveSurface () { return &ss; }

	//! return integration rule class
	bool UseNodalIntegration() { return true; }

	//! build the matrix profile for use in the stiffness matrix
	void BuildMatrixProfile(FEGlobalMatrix& K);

public:
	FEStickySurface	ss;	//!< slave surface
	FEStickySurface	ms;	//!< master surface

	int nse;	//!< number of slave elements
	int nme;	//!< number of master elements

	double		m_atol;		//!< augmentation tolerance
	double		m_eps;		//!< penalty scale factor
	double		m_stol;		//!< search tolerance
	int			m_naugmax;	//!< maximum nr of augmentations
	int			m_naugmin;	//!< minimum nr of augmentations
	double		m_tmax;		//!< max traction
	double		m_snap;		//!< snap tolerance

	DECLARE_PARAMETER_LIST();
};
