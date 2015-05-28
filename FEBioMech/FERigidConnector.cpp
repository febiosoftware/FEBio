//
//  FERigidConnector.cpp
//  FEBioMech
//
//  Created by Gerard Ateshian on 5/27/15.
//  Copyright (c) 2015 febio.org. All rights reserved.
//

#include "FERigidConnector.h"

//-----------------------------------------------------------------------------
BEGIN_PARAMETER_LIST(FERigidConnector, FENLConstraint);
ADD_PARAMETER(m_nRBa, FE_PARAM_INT   , "body_a"        );
ADD_PARAMETER(m_nRBb, FE_PARAM_INT   , "body_b"        );
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FERigidConnector::FERigidConnector(FEModel* pfem) : FENLConstraint(pfem) {
    m_F = m_M = vec3d(0, 0, 0);
};

//-----------------------------------------------------------------------------
FERigidConnector::~FERigidConnector() {}