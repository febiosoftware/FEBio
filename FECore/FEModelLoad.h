#pragma once
#include "FEModelComponent.h"
#include "FEGlobalVector.h"
#include "FETypes.h"

//-----------------------------------------------------------------------------
class FESolver;

//-----------------------------------------------------------------------------
//! This class is the base class for all classes that affect the state of the model
//! and contribute directly to the residual and the global stiffness matrix. This
//! includes most boundary loads, body loads, contact, etc.
class FEModelLoad : public FEModelComponent
{
public:
	//! constructor
	FEModelLoad(SUPER_CLASS_ID sid, FEModel* pfem);

public:
	// all classes derived from this base class must implement
	// the following functions.

	//! evaluate the contribution to the residual
	virtual void Residual(FEGlobalVector& R, const FETimePoint& tp);

	//! evaluate the contribution to the global stiffness matrix
	virtual void StiffnessMatrix(FESolver* psolver, const FETimePoint& tp);
};

//-----------------------------------------------------------------------------
//! an axial force between two rigid bodies
class FERigidAxialForce : public FEModelLoad
{
public:
	//! constructor
	FERigidAxialForce(FEModel* pfem);

	//! initialization
	bool Init();

	//! Serialization
	void Serialize(DumpFile& ar);

	//! Residual
	void Residual(FEGlobalVector& R, const FETimePoint& tp);

	//! Stiffness matrix
	void StiffnessMatrix(FESolver* psolver, const FETimePoint& tp);

public:
	int		m_ida, m_idb;		//!< rigid body ID's
	vec3d	m_ra0, m_rb0;		//!< coordinates of attachements in reference state
	double	m_s;				//!< scale factor
	bool	m_brelative;		//!< if active, the ra0 and rb0 are relative w.r.t. the COM

	DECLARE_PARAMETER_LIST();
};


//-----------------------------------------------------------------------------
//! rigid body force

class FERigidBodyForce : public FEModelLoad
{
public:
	FERigidBodyForce(FEModel* pfem);

	//! get the current force value
	double Value();

	//! Serialization
	void Serialize(DumpFile& ar);

	//! Residual
	void Residual(FEGlobalVector& R, const FETimePoint& tp);

	//! Stiffness matrix
	void StiffnessMatrix(FESolver* psolver, const FETimePoint& tp);

public:
	int		ntype;	//!< type of force (0=loadcurve, 1=target)
	int		id;		// rigid body id
	int		bc;		// force direction
	int		lc;		// load curve number
	double	sf;		// scale factor
	bool	m_bfollow;	//!< follower force if true
};
