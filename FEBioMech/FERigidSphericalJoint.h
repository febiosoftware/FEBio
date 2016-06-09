// FERigidSphericalJoint.h: interface for the FERigidSphericalJoint class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FERigidSphericalJoint_H__22EBB207_A1C7_465E_B985_A2140CBA9BDB__INCLUDED_)
#define AFX_FERigidSphericalJoint_H__22EBB207_A1C7_465E_B985_A2140CBA9BDB__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "FECore/vec3d.h"
#include "FERigidConnector.h"

//-----------------------------------------------------------------------------
//! The FERigidSphericalJoint class implements a rigid spherical joint that
//! functions under static and dynamic conditions. This joint allows the
//! user to connect two rigid bodies at a point in space.

class FERigidSphericalJoint : public FERigidConnector
{
public:
    //! constructor
    FERigidSphericalJoint(FEModel* pfem);
    
    //! destructor
    virtual ~FERigidSphericalJoint() {}
    
    //! initialization
    bool Init();
    
    //! calculates the joint forces
    void Residual(FEGlobalVector& R, const FETimeInfo& tp);
    
    //! calculates the joint stiffness
    void StiffnessMatrix(FESolver* psolver, const FETimeInfo& tp);
    
    //! calculate Lagrangian augmentation
    bool Augment(int naug, const FETimeInfo& tp);
    
    //! serialize data to archive
    void Serialize(DumpStream& ar);
    
    //! update state
    void Update(const FETimeInfo& tp);
    
    //! Reset data
    void Reset();
    
public: // parameters
    double	m_atol;     //! augmented Lagrangian tolerance
    double  m_gtol;     //! augmented Lagrangian gap tolerance
    double  m_qtol;     //! augmented Lagrangian angular gap tolerance
    double	m_eps;      //! penalty factor for constraining force
    double	m_ups;      //! penalty factor for constraining moment
    vec3d	m_q0;       //! initial position of joint
    int     m_naugmin;  //! minimum number of augmentations
    int     m_naugmax;  //! maximum number of augmentations
    bool    m_bq;           //! flag for prescribing rotation
    double  m_qpx;          //! prescribed rotation x-component
    double  m_qpy;          //! prescribed rotation y-component
    double  m_qpz;          //! prescribed rotation z-component
    double  m_Mpx;          //! prescribed moment along x
    double  m_Mpy;          //! prescribed moment along y
    double  m_Mpz;          //! prescribed moment along z

protected:
    vec3d	m_qa0;      //! initial relative position vector of joint w.r.t. A
    vec3d	m_qb0;      //! initial relative position vector of joint w.r.t. B
   
    vec3d	m_e0[3];        //! initial joint basis
    vec3d	m_ea0[3];       //! initial joint basis w.r.t. A
    vec3d	m_eb0[3];       //! initial joint basis w.r.t. B
    
    vec3d	m_L;        //! Lagrange multiplier for constraining force
    vec3d	m_U;        //! Lagrange multiplier for constraining moment
    
protected:
    bool	m_binit;
    
    DECLARE_PARAMETER_LIST();
};

#endif // !defined(AFX_FERigidSphericalJoint_H__22EBB207_A1C7_465E_B985_A2140CBA9BDB__INCLUDED_)
