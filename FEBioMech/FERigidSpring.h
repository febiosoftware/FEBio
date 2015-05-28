//
//  FERigidSpring.h
//  FEBioMech
//
//  Created by Gerard Ateshian on 5/11/15.
//  Copyright (c) 2015 febio.org. All rights reserved.
//

#ifndef __FEBioMech__FERigidSpring__
#define __FEBioMech__FERigidSpring__

#include "FECore/vec3d.h"
#include "FECore/DumpFile.h"
#include "FERigidConnector.h"

//-----------------------------------------------------------------------------
//! The FERigidSpring class implements a linear spring that connects
//! two rigid bodies at arbitrary points (not necessarily nodes).

class FERigidSpring : public FERigidConnector
{
public:
    //! constructor
    FERigidSpring(FEModel* pfem);
    
    //! destructor
    ~FERigidSpring() {}
    
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
    vec3d	m_a0;       //! initial absolute position vector of spring on body A
    vec3d	m_b0;       //! initial absolute position vector of spring on body B
    vec3d	m_qa0;      //! initial relative position vector of spring on body A
    vec3d	m_qb0;      //! initial relative position vector of spring on body B
    
    double	m_k;        //! spring constant
    double  m_L0;       //! spring resting length
    
protected:
    int		m_nID;      //!< ID of rigid joint
    bool	m_binit;
    double  m_alpha;    //! alpha from solver
    
    DECLARE_PARAMETER_LIST();
};

#endif /* defined(__FEBioMech__FERigidSpring__) */
