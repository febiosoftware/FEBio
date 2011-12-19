// FERigidJoint.h: interface for the FERigidJoint class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FERIGIDJOINT_H__22EBB207_A1C7_465E_B985_A2140CBA9BDB__INCLUDED_)
#define AFX_FERIGIDJOINT_H__22EBB207_A1C7_465E_B985_A2140CBA9BDB__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "FECore/vec3d.h"
#include "NumCore/vector.h"
#include "FECore/DumpFile.h"
#include "FECore/FEParameterList.h"

class FEM;

//-----------------------------------------------------------------------------
//! The FERigidJoint class implements a rigid joint. The rigid joint allows the
//! user to connect two rigid bodies at a point in space.

class FERigidJoint : public FEParamContainer
{
public:
	//! constructor
	FERigidJoint(FEM* pfem);

	//! destructor
	virtual ~FERigidJoint();

	//! create a shallow copy
	void ShallowCopy(FERigidJoint& rj)
	{
		m_nRBa = rj.m_nRBa;
		m_nRBb = rj.m_nRBb;

		m_q0  = rj.m_q0;
		m_qa0 = rj.m_qa0;
		m_qb0 = rj.m_qb0;

		m_F = rj.m_F;
		m_L = rj.m_L;

		m_eps  = rj.m_eps;
		m_atol = rj.m_atol;
	}

	//! calculates the joint forces
	void JointForces(vector<double>& R);

	//! calculates the joint stiffness
	void JointStiffness();

	//! calculate Lagrangian augmentation
	bool Augment();

	//! serialize data to archive
	void Serialize(DumpFile& ar);

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
	FEM*	m_pfem;	//!< FEM class to which this rigid joint belongs

	int	m_nID;	//!< ID of rigid joint

	DECLARE_PARAMETER_LIST();
};

#endif // !defined(AFX_FERIGIDJOINT_H__22EBB207_A1C7_465E_B985_A2140CBA9BDB__INCLUDED_)
