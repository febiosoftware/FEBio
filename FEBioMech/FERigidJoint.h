#pragma once
#include "FECore/vec3d.h"
#include "FERigidConnector.h"

//-----------------------------------------------------------------------------
//! The FERigidJoint class implements a rigid joint. The rigid joint allows the
//! user to connect two rigid bodies at a point in space.

class FERigidJoint : public FERigidConnector
{
public:
	//! constructor
	FERigidJoint(FEModel* pfem);

	//! destructor
	virtual ~FERigidJoint() {}

	//! initialization
	bool Init() override;

	//! calculates the joint forces
	void Residual(FEGlobalVector& R, const FETimeInfo& tp) override;

	//! calculates the joint stiffness
	void StiffnessMatrix(FESolver* psolver, const FETimeInfo& tp) override;

	//! calculate Lagrangian augmentation
	bool Augment(int naug, const FETimeInfo& tp) override;

	//! serialize data to archive
	void Serialize(DumpStream& ar) override;

	//! update state
	void Update() override;

	//! Reset data
	void Reset() override;

public:
	vec3d	m_q0;		//! initial position of joint
	vec3d	m_qa0;	//! initial relative position vector of joint w.r.t. A
	vec3d	m_qb0;	//! initial relative position vector of joint w.r.t. B

	vec3d	m_F;		//!< constraining force
	vec3d	m_L;		//!< Lagrange multiplier
	double	m_eps;		//!< penalty factor
	double	m_atol;		//!< augmented Lagrangian tolerance
	bool	m_blaugon;	//!< augmented Lagrangian flag

protected:
	int		m_nID;	//!< ID of rigid joint
	bool	m_binit;

	DECLARE_FECORE_CLASS();
};
