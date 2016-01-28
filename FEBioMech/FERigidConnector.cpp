//
//  FERigidConnector.cpp
//  FEBioMech
//
//  Created by Gerard Ateshian on 5/27/15.
//  Copyright (c) 2015 febio.org. All rights reserved.
//

#include "FERigidConnector.h"
#include "FECore/FEModel.h"

//-----------------------------------------------------------------------------
BEGIN_PARAMETER_LIST(FERigidConnector, FENLConstraint);
ADD_PARAMETER(m_nRBa, FE_PARAM_INT   , "body_a"        );
ADD_PARAMETER(m_nRBb, FE_PARAM_INT   , "body_b"        );
END_PARAMETER_LIST();

int FERigidConnector::m_ncount = 0;

//-----------------------------------------------------------------------------
FERigidConnector::FERigidConnector(FEModel* pfem) : FENLConstraint(pfem) {
    m_F = m_M = vec3d(0, 0, 0);
};

//-----------------------------------------------------------------------------
FERigidConnector::~FERigidConnector() {}

//-----------------------------------------------------------------------------
void FERigidConnector::BuildMatrixProfile(FEGlobalMatrix& M)
{
	FERigidSystem& rigid = *GetFEModel()->GetRigidSystem();
    vector<int> lm(12);
                    
    int* lm1 = rigid.Object(m_nRBa)->m_LM;
    int* lm2 = rigid.Object(m_nRBb)->m_LM;
                    
    for (int j=0; j<6; ++j) lm[j  ] = lm1[j];
    for (int j=0; j<6; ++j) lm[j+6] = lm2[j];
    M.build_add(lm);
}

//-----------------------------------------------------------------------------
void FERigidConnector::Serialize(DumpStream& ar)
{
	FENLConstraint::Serialize(ar);
	if (ar.IsSaving())
	{
        ar << m_nID;
        ar << m_nRBa << m_nRBb;
        ar << m_F << m_M;
	}
	else
	{
        ar >> m_nID;
        ar >> m_nRBa >> m_nRBb;
        ar >> m_F >> m_M;
	}
}
