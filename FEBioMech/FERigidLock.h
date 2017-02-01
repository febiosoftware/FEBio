//
//  FERigidLock.hpp
//  FEBioMech
//
//  Created by Gerard Ateshian on 1/31/17.
//  Copyright © 2017 febio.org. All rights reserved.
//

#ifndef FERigidLock_hpp
#define FERigidLock_hpp

#include "FECore/vec3d.h"
#include "FERigidConnector.h"

//-----------------------------------------------------------------------------
//! The FERigidLock class implements a locking joint. This rigid joint
//! allows the user to connect two rigid bodies at a point in space
//! and prevents any relative motion.

class FERigidLock : public FERigidConnector
{
public:
    //! constructor
    FERigidLock(FEModel* pfem);
    
    //! destructor
    virtual ~FERigidLock() {}
    
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
    double	m_atol;	//! augmented Lagrangian tolerance
    double  m_gtol; //! augmented Lagrangian gap tolerance
    double  m_qtol; //! augmented Lagrangian angular gap tolerance
    int     m_naugmin;  //! minimum number of augmentations
    int     m_naugmax;  //! maximum number of augmentations
    vec3d	m_q0;	//! initial position of joint
    double	m_eps;	//! penalty factor for constraining force
    double	m_ups;	//! penalty factor for constraining moment
    
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

#endif /* FERigidLock_hpp */
