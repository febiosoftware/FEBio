// FERigidJoint.h: interface for the FERigidJoint class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FERIGIDJOINT_H__22EBB207_A1C7_465E_B985_A2140CBA9BDB__INCLUDED_)
#define AFX_FERIGIDJOINT_H__22EBB207_A1C7_465E_B985_A2140CBA9BDB__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "FECore/vec3d.h"
#include "FECore/DumpFile.h"
#include "FECore/FENLConstraint.h"

//-----------------------------------------------------------------------------
//! The FERigidJoint class implements a rigid joint. The rigid joint allows the
//! user to connect two rigid bodies at a point in space.

class FERigidJoint : public FENLConstraint
{
public:
	//! constructor
	FERigidJoint(FEModel* pfem);

	//! destructor
	virtual ~FERigidJoint() {}

	//! create a shallow copy
	void ShallowCopy(FERigidJoint& rj);

public:

	//! initialization \todo Find a use for this
	void Init() {}

	//! calculates the joint forces
	void Residual(FEGlobalVector& R);

	//! calculates the joint stiffness
	void StiffnessMatrix(FESolver* psolver);

	//! calculate Lagrangian augmentation
	bool Augment(int naug);

	//! serialize data to archive
	void Serialize(DumpFile& ar);

	//! update state
	void Update();

public:
	int	m_nRBa;		//!< rigid body A that the joint connects
	int	m_nRBb;		//!< rigid body B that the joint connects

	vec3d	m_q0;		//! initial position of joint
	vec3d	m_qa0;	//! initial relative position vector of joint w.r.t. A
	vec3d	m_qb0;	//! initial relative position vector of joint w.r.t. B

	vec3d	m_F;		//! constraining force
	vec3d	m_L;		//! Lagrange multiplier
	double	m_eps;	//! penalty factor
	double	m_atol;	//! augmented Lagrangian tolerance

protected:
	int	m_nID;	//!< ID of rigid joint
	DECLARE_PARAMETER_LIST();
};

#endif // !defined(AFX_FERIGIDJOINT_H__22EBB207_A1C7_465E_B985_A2140CBA9BDB__INCLUDED_)
