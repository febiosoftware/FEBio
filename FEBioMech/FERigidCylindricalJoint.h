//
//  FERigidCylindricalJoint.h
//  FEBioMech
//
//  Created by Gerard Ateshian on 5/9/15.
//  Copyright (c) 2015 febio.org. All rights reserved.
//

#ifndef __FEBioMech__FERigidCylindricalJoint__
#define __FEBioMech__FERigidCylindricalJoint__

#include "FECore/vec3d.h"
#include "FECore/DumpFile.h"
#include "FECore/LoadCurve.h"
#include "FERigidConnector.h"

//-----------------------------------------------------------------------------
//! The FERigidCylindricalJoint class implements a cylindrical joint. The rigid joint
//! allows the user to connect two rigid bodies at a point in space
//! and allow rotation about, and translation along, a single prescribed axis.

class FERigidCylindricalJoint : public FERigidConnector
{
public:
    //! constructor
    FERigidCylindricalJoint(FEModel* pfem);
    
    //! destructor
    virtual ~FERigidCylindricalJoint() {}
    
    //! initialization
    bool Init();
    
    //! calculates the joint forces
    void Residual(FEGlobalVector& R, const FETimePoint& tp);
    
    //! calculates the joint stiffness
    void StiffnessMatrix(FESolver* psolver, const FETimePoint& tp);
    
    //! calculate Lagrangian augmentation
    bool Augment(int naug, const FETimePoint& tp);
    
    //! serialize data to archive
    void Serialize(DumpFile& ar);
    
    //! create a shallow copy
    void ShallowCopy(DumpStream& dmp, bool bsave);
    
    //! update state
    void Update(const FETimePoint& tp);
    
    //! Reset data
    void Reset();
    
public: // parameters
    double	m_atol;	//! augmented Lagrangian tolerance
    double  m_gtol; //! augmented Lagrangian gap tolerance
    double  m_qtol; //! augmented Lagrangian angular gap tolerance
    int     m_naugmin;  //! minimum number of augmentations
    int     m_naugmax;  //! maximum number of augmentations
    double	m_eps;	//! penalty factor for constraining force
    double	m_ups;	//! penalty factor for constraining moment
    vec3d	m_q0;	//! initial position of joint
    double  m_dp;   //! prescribed translation
    double  m_qp;   //! prescribed rotation
    bool    m_bd;   //! flag for prescribing translation
    bool    m_bq;   //! flag for prescribing rotation
    double  m_Fp;   //! prescribed force
    double  m_Mp;   //! prescribed moment

protected:
    vec3d	m_qa0;	//! initial relative position vector of joint w.r.t. A
    vec3d	m_qb0;	//! initial relative position vector of joint w.r.t. B
    
    vec3d	m_e0[3];	//! initial joint basis
    vec3d	m_ea0[3];	//! initial joint basis w.r.t. A
    vec3d	m_eb0[3];	//! initial joint basis w.r.t. B
    
    vec3d	m_L;	//! Lagrange multiplier for constraining force
    vec3d	m_U;	//! Lagrange multiplier for constraining moment
    
protected:
    bool	m_binit;

	DECLARE_PARAMETER_LIST();
};

#endif /* defined(__FEBioMech__FERigidCylindricalJoint__) */
