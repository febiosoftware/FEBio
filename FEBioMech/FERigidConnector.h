//
//  FERigidConnector.h
//  FEBioMech
//
//  Created by Gerard Ateshian on 5/27/15.
//  Copyright (c) 2015 febio.org. All rights reserved.
//

#ifndef __FEBioMech__FERigidConnector__
#define __FEBioMech__FERigidConnector__

#include "FECore/FENLConstraint.h"

//-----------------------------------------------------------------------------
//! This is a virtual class for all rigid connectors, including
//! spherical, revolute, prismatic and cylindrical joints, as well
//! as springs and dampers that connect rigid bodies.

class FERigidConnector : public FENLConstraint
{
public:
    //! constructor
    FERigidConnector(FEModel* pfem);
    
    //! destructor
    virtual ~FERigidConnector();
    
    int     GetConnectorID() { return m_nID; }

	//! build connectivity for matrix profile
	void BuildMatrixProfile(FEGlobalMatrix& M);

	//! serialization
	void Serialize(DumpStream& ar);
    
public:
    int	m_nRBa;		//!< rigid body A that the connector connects
    int	m_nRBb;		//!< rigid body B that the connector connects
    
    vec3d	m_F;	//! constraining force
    vec3d	m_M;	//! constraining moment
    
protected:
    int		m_nID;	//!< ID of rigid connector
    
    static int	m_ncount;	//!< used to create unique ID's for the nonlinear constraints
    
    DECLARE_PARAMETER_LIST();
};

#endif /* defined(__FEBioMech__FERigidConnector__) */
