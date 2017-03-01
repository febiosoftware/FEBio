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
#include "FEObject.h"

//-----------------------------------------------------------------------------
class FEModel;
class FERigidBodyDisplacement;

//-----------------------------------------------------------------------------
//! rigid body class
//! \todo perhaps the rigid body should store a list of domains it uses.
//!       That way, we can have multiple domains per RB using multiple 
//!       materials.
class FERigidBody : public FEObject 
{
public:
	// Constructor
	FERigidBody(FEModel* pfem);

	//! desctructor
	virtual ~FERigidBody();

	//! Set the center of mass directly
	void SetCOM(vec3d rc);

	//! Update total mass and center of mass
	void UpdateCOM();

	//! reset rigid body data
	void Reset();

	//! initialize data
	void Init();

	//! serialize data to archive
	void Serialize(DumpStream& ar);

	//! get the material ID
	int GetMaterialID() { return m_mat; }
    
    //! incremental compound rotation from Cayley transform
    vec3d CayleyIncrementalCompoundRotation();

public:
	const quatd& GetRotation() const { return m_qt; }
	quatd GetPreviousRotation() const { return m_qp; }
	void SetRotation(const quatd& q)
	{
		m_qt = q;
		m_qt.GetEuler(m_euler.x, m_euler.y, m_euler.z);
	}

public:
	int		m_nID;		//!< ID of rigid body
	int		m_mat;		//!< material ID (TODO: Since rigid bodies can have multiple materials, I want to remove this)
	double	m_mass;		//!< total mass of rigid body
    mat3ds  m_moi;      //!< mass moment of inertia about center of mass
	vec3d	m_Fr, m_Mr;	//!< reaction force and torque
	vec3d	m_Fp, m_Mp;	//!< reaction force and torque at the end of the previous step

	vec3d	m_r0;	//!< initial position of rigid body
	vec3d	m_rp;	//!< previous position of rigid body
	vec3d	m_rt;	//!< current position of rigid body
    
	vec3d	m_vp;	//!< previous velocity of rigid body
	vec3d	m_vt;	//!< current velocity of rigid body
    
	vec3d	m_ap;	//!< previous acceleration of rigid body
	vec3d	m_at;	//!< current acceleration of rigid body

	quatd	m_qp;	//!< previous orientation of rigid body

private:
	// TODO: This is a hack!I only need this so I can access the euler angles directly from
	//       the optimization module. Need to figure out a better way.
	quatd	m_qt;		//!< current orientation of rigid body
	vec3d	m_euler;	//!< Euler angles of rotation 

public:
    vec3d   m_wp;   //!< previous angular velocity of rigid body
    vec3d   m_wt;   //!< current angular velocity of rigid body
    
    vec3d   m_alp;  //!< previous angular acceleration of rigid body
    vec3d   m_alt;  //!< current angular acceleration of rigid body
    
	int		m_BC[6];	//!< DOF types
	int		m_LM[6];	//!< dof equation numbers
	double	m_Up[6];	//!< previous displacement/rotation vector
	double	m_Ut[6];	//!< total displacement/rotation vector
	double	m_du[6];	//!< incremental displacement vector
	double	m_dul[6];	//!< displacement in local coordinates system

    bool    m_bpofr;    //!< flag for all or none of rotation dofs prescribed/fixed
    
public:
	FERigidBodyDisplacement*	m_pDC[6];	//!< active displacement constraints
	FERigidBody*	m_prb;	//!< parent rigid body

public:
	DECLARE_PARAMETER_LIST();
};

#endif // !defined(AFX_FERIGIDBODY_H__2C1FB6E7_60F0_46E9_94D0_3B4D07EAC5CF__INCLUDED_)
