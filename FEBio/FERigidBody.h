// FERigidBody.h: interface for the FERigidBody class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FERIGIDBODY_H__2C1FB6E7_60F0_46E9_94D0_3B4D07EAC5CF__INCLUDED_)
#define AFX_FERIGIDBODY_H__2C1FB6E7_60F0_46E9_94D0_3B4D07EAC5CF__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "FE_enum.h"
#include "FEBioLib/vec3d.h"
#include "FEBioLib/quatd.h"
#include "FEBoundaryCondition.h"
#include "DumpFile.h"

class FEM;


//-----------------------------------------------------------------------------
//! rigid body force

class FERigidBodyForce : public FEBoundaryCondition
{
public:
	int		id;	// rigid body id
	int		bc;	// force direction
	int		lc;	// load curve number
	double	sf;	// scale factor
};

//-----------------------------------------------------------------------------
//! rigid body displacement

class FERigidBodyDisplacement : public FEBoundaryCondition
{
public:
	virtual ~FERigidBodyDisplacement(){}

public:
	int		id;	// rigid body id
	int		bc;	// displacement direction
	int		lc;	// load curve number
	double	sf;	// scale factor
};

//-----------------------------------------------------------------------------
//! rigid body fixed constraint
//! This constraint is used to fix a rigid body

class FERigidFixedConstraint : public FEBoundaryCondition
{
public:
	int		id;	// rigid body id
	int		bc;	// displacement direction
};

//-----------------------------------------------------------------------------
//! rigid body class

class FERigidBody  
{
public:
	FERigidBody();
	virtual ~FERigidBody();

	void AttachToFEM(FEM* pfem) { m_pfem = pfem; }

	void Update();

	//! serialize data to archive
	void Serialize(DumpFile& ar);

public:
	int		m_nID;		// ID of rigid body
	int		m_mat;		// material ID
	double	m_mass;		// total mass of rigid body
	vec3d	m_Fr, m_Mr;	// reaction force and torque

	vec3d	m_r0;	// initial position of rigid body
	vec3d	m_rp;	// previous position of rigid body
	vec3d	m_rt;	// current position of rigid body

	quatd	m_qp;	// previous orientation of rigid body
	quatd	m_qt;	// current orientation of rigid body

	int		m_LM[6];	// equation numbers
	double	m_Up[6];	// previous displacement/rotation vector
	double	m_Ut[6];	// total displacement/rotation vector
	double	m_du[6];	// incremental displacement vector
	double	m_dul[6];	// displacement in local coordinates system

	FERigidBodyDisplacement*	m_pDC[6];	// active displacement constraints

	FERigidBody*	m_prb;	//!< parent rigid body

	bool	m_bActive;	// activation flag

private:
	FEM*	m_pfem;
};

#endif // !defined(AFX_FERIGIDBODY_H__2C1FB6E7_60F0_46E9_94D0_3B4D07EAC5CF__INCLUDED_)
