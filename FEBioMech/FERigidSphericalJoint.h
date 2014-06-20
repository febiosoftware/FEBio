// FERigidSphericalJoint.h: interface for the FERigidSphericalJoint class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FERigidSphericalJoint_H__22EBB207_A1C7_465E_B985_A2140CBA9BDB__INCLUDED_)
#define AFX_FERigidSphericalJoint_H__22EBB207_A1C7_465E_B985_A2140CBA9BDB__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "FECore/vec3d.h"
#include "FECore/DumpFile.h"
#include "FECore/FENLConstraint.h"

//-----------------------------------------------------------------------------
//! The FERigidSphericalJoint class implements a rigid spherical joint that
//! functions under static and dynamic conditions. This joint allows the
//! user to connect two rigid bodies at a point in space.

class FERigidSphericalJoint : public FENLConstraint
{
public:
	//! constructor
	FERigidSphericalJoint(FEModel* pfem);
    
	//! destructor
	virtual ~FERigidSphericalJoint() {}
    
	//! initialization
	bool Init();
    
	//! calculates the joint forces
	void Residual(FEGlobalVector& R);
    
	//! calculates the joint stiffness
	void StiffnessMatrix(FESolver* psolver);
    
	//! calculate Lagrangian augmentation
	bool Augment(int naug);
    
	//! serialize data to archive
	void Serialize(DumpFile& ar);
    
	//! create a shallow copy
	void ShallowCopy(DumpStream& dmp, bool bsave);
    
	//! update state
	void Update();
    
	//! Reset data
	void Reset();
    
public:
	int	m_nRBa;		//!< rigid body A that the joint connects
	int	m_nRBb;		//!< rigid body B that the joint connects
    
	vec3d	m_q0;		//! initial position of joint
	vec3d	m_qa0;	//! initial relative position vector of joint w.r.t. A
	vec3d	m_qb0;	//! initial relative position vector of joint w.r.t. B
    
	vec3d	m_F;		//! constraining force
	vec3d	m_L;		//! Lagrange multiplier
	double	m_eps;      //! penalty factor
	double	m_atol;     //! augmented Lagrangian tolerance
    double  m_gtol;     //! augmented Lagrangian gap tolerance
    int     m_naugmin;  //! minimum number of augmentations
    int     m_naugmax;  //! maximum number of augmentations
    
protected:
	int		m_nID;	//!< ID of rigid joint
	bool	m_binit;
    double  m_alpha;    //! alpha from solver
    
	DECLARE_PARAMETER_LIST();
};

#endif // !defined(AFX_FERigidSphericalJoint_H__22EBB207_A1C7_465E_B985_A2140CBA9BDB__INCLUDED_)
