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
#include "BC.h"
#include "DumpFile.h"
#include "FEModel.h"
#include "FEObject.h"

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

	//! update solution
	void Update(std::vector<double>& Ui, std::vector<double>& ui);

	//! serialize data to archive
	void Serialize(DumpFile& ar);

	//! shallow copy
	void ShallowCopy(DumpStream& dmp, bool bsave);

	//! State managment
	bool IsActive() { return m_bActive; }
	void Activate(bool b) { m_bActive = b; }

	//! get the material ID
	int GetMaterialID() { return m_mat; }

public:
	int		m_nID;		//!< ID of rigid body
	int		m_mat;		//!< material ID
	double	m_mass;		//!< total mass of rigid body
    mat3ds  m_moi;      //!< mass moment of inertia about center of mass
	vec3d	m_Fr, m_Mr;	//!< reaction force and torque

	vec3d	m_r0;	//!< initial position of rigid body
	vec3d	m_rp;	//!< previous position of rigid body
	vec3d	m_rt;	//!< current position of rigid body

	quatd	m_qp;	//!< previous orientation of rigid body
	quatd	m_qt;	//!< current orientation of rigid body

	vec3d	m_vt;	//!< linear velocity
	vec3d	m_wt;	//!< angular velocity

	int		m_BC[6];	//!< dof constrains (0=free, -1=fixed, >0 = prescribed)
	int		m_LM[6];	//!< dof equation numbers
	double	m_Up[6];	//!< previous displacement/rotation vector
	double	m_Ut[6];	//!< total displacement/rotation vector
	double	m_du[6];	//!< incremental displacement vector
	double	m_dul[6];	//!< displacement in local coordinates system

public:
	FERigidBodyDisplacement*	m_pDC[6];	//!< active displacement constraints
	FERigidBody*	m_prb;	//!< parent rigid body

protected:
	bool	m_bActive;	// activation flag

public:
	DECLARE_PARAMETER_LIST();
};

#endif // !defined(AFX_FERIGIDBODY_H__2C1FB6E7_60F0_46E9_94D0_3B4D07EAC5CF__INCLUDED_)
