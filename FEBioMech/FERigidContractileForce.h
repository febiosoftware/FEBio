//
//  FERigidContractileForce.h
//  FEBioMech
//
//  Created by Gerard Ateshian on 5/23/15.
//  Copyright (c) 2015 febio.org. All rights reserved.
//

#ifndef __FEBioMech__FERigidContractileForce__
#define __FEBioMech__FERigidContractileForce__

#include "FECore/vec3d.h"
#include "FECore/DumpFile.h"
#include "FERigidConnector.h"

//-----------------------------------------------------------------------------
//! The FERigidContractileForce class implements a contractile force between
//! arbitrary points (not necessarily nodes) on two rigid bodies.

class FERigidContractileForce : public FERigidConnector
{
public:
    //! constructor
    FERigidContractileForce(FEModel* pfem);
    
    //! destructor
    ~FERigidContractileForce() {}
    
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
    vec3d	m_a0;       //! initial absolute position vector of insertion on body A
    vec3d	m_b0;       //! initial absolute position vector of insertion on body B
    vec3d	m_qa0;      //! initial relative position vector of insertion on body A
    vec3d	m_qb0;      //! initial relative position vector of insertion on body B
    
    double	m_f0;       //! contractile force
    
protected:
    bool	m_binit;
    
    DECLARE_PARAMETER_LIST();
};

#endif /* defined(__FEBioMech__FERigidContractileForce__) */
