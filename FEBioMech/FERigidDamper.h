//
//  FERigidDamper.h
//  FEBioMech
//
//  Created by Gerard Ateshian on 5/12/15.
//  Copyright (c) 2015 febio.org. All rights reserved.
//

#ifndef __FEBioMech__FERigidDamper__
#define __FEBioMech__FERigidDamper__

#include "FECore/vec3d.h"
#include "FECore/DumpFile.h"
#include "FERigidConnector.h"

//-----------------------------------------------------------------------------
//! The FERigidDamper class implements a linear damper that connects
//! two rigid bodies at arbitrary points (not necessarily nodes).

class FERigidDamper : public FERigidConnector
{
public:
    //! constructor
    FERigidDamper(FEModel* pfem);
    
    //! destructor
    ~FERigidDamper() {}
    
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
    
public:
    vec3d	m_a0;       //! initial absolute position vector of spring on body A
    vec3d	m_b0;       //! initial absolute position vector of spring on body B
    vec3d	m_qa0;      //! initial relative position vector of spring on body A
    vec3d	m_qb0;      //! initial relative position vector of spring on body B
    
    double	m_c;        //! damping constant
    
protected:
    int		m_nID;      //!< ID of rigid joint
    bool	m_binit;
    
    DECLARE_PARAMETER_LIST();
};

#endif /* defined(__FEBioMech__FERigidDamper__) */
