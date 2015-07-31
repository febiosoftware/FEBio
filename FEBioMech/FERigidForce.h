#pragma once
#include "FECore/FEModelLoad.h"

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
//! TODO: I'd like to split this class into two classes: one handling the case
//!       were the force is const, and one where the force is a follower force.
//!       Perhaps I can derive the const force from FENodalLoad since it applies
//!       a force directly to the rigid "node".
class FERigidBodyForce : public FEModelLoad
{
public:
	enum { RAMP, TARGET };	// values for m_ntype

public:
	FERigidBodyForce(FEModel* pfem);

	//! initialization
	bool Init();

	//! get the current force value
	double Value();

	//! Serialization
	void Serialize(DumpFile& ar);

	//! Residual
	void Residual(FEGlobalVector& R, const FETimePoint& tp);

	//! Stiffness matrix
	void StiffnessMatrix(FESolver* psolver, const FETimePoint& tp);

public:
	int		m_ntype;	//!< type of force (0=loadcurve, 1=target)
	int		id;		// rigid body id
	int		bc;		// force direction
	int		lc;		// load curve number
	double	sf;		// scale factor
	bool	m_bfollow;	//!< follower force if true
};
