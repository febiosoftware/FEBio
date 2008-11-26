// FERigidBody.h: interface for the FERigidBody class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FERIGIDBODY_H__2C1FB6E7_60F0_46E9_94D0_3B4D07EAC5CF__INCLUDED_)
#define AFX_FERIGIDBODY_H__2C1FB6E7_60F0_46E9_94D0_3B4D07EAC5CF__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "FE_enum.h"
#include "vec3d.h"
#include "quatd.h"

class FEM;

class FERigidBody  
{
public:
	FERigidBody();
	virtual ~FERigidBody();

	void AttachToFEM(FEM* pfem) { m_pfem = pfem; }

	void Update();

public:
	int		m_nID;	// ID of rigid body
	int		m_mat;		// material ID
	double	m_mass;		// total mass of rigid body
	vec3d	m_Fr, m_Mr;	// reaction force and torque

	vec3d	m_r0;	// initial position of rigid body
	vec3d	m_rp;	// previous position of rigid body
	vec3d	m_rt;	// current position of rigid body

	quatd	m_qp;	// previous orientation of rigid body
	quatd	m_qt;	// current orientation of rigid body

	int		m_bc[6];	// constraints
	int		m_LM[6];	// equation numbers
	double	m_Up[6];	// previous displacement/rotation vector
	double	m_Ut[6];	// total displacement/rotation vector
	double	m_du[6];	// incremental displacement vector

private:
	FEM*	m_pfem;
};

#endif // !defined(AFX_FERIGIDBODY_H__2C1FB6E7_60F0_46E9_94D0_3B4D07EAC5CF__INCLUDED_)
