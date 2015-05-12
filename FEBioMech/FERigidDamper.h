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
#include "FECore/FENLConstraint.h"

//-----------------------------------------------------------------------------
//! The FERigidDamper class implements a linear damper that connects
//! two rigid bodies at arbitrary points (not necessarily nodes).

class FERigidDamper : public FENLConstraint
{
public:
    //! constructor
    FERigidDamper(FEModel* pfem);
    
    //! destructor
    ~FERigidDamper() {}
    
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
    int	m_nRBa;         //!< rigid body A that the spring connects
    int	m_nRBb;         //!< rigid body B that the spring connects
    
    vec3d	m_a0;       //! initial absolute position vector of spring on body A
    vec3d	m_b0;       //! initial absolute position vector of spring on body B
    vec3d	m_qa0;      //! initial relative position vector of spring on body A
    vec3d	m_qb0;      //! initial relative position vector of spring on body B
    
    vec3d	m_F;		//! constraining force
    double	m_c;        //! damping constant
    
protected:
    int		m_nID;      //!< ID of rigid joint
    bool	m_binit;
    double  m_alpha;    //! alpha from solver
    double  m_beta;     //! beta from solver
    double  m_gamma;    //! gamma from solver
    
    DECLARE_PARAMETER_LIST();
};

#endif /* defined(__FEBioMech__FERigidDamper__) */
