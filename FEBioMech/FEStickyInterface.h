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
    void GetContactTraction(int nface, vec3d& pt);
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
	bool Init() override;

	//! interface activation
	void Activate() override;

	//! projects slave nodes onto master nodes
	void ProjectSurface(FEStickySurface& ss, FEStickySurface& ms, bool bmove = false);

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
	void Update(int niter, const FETimeInfo& tp) override;

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

	DECLARE_FECORE_CLASS();
};
