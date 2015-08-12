//
//  FEContinuousFiberDistribution.cpp
//  FEBioXCode4
//
//  Created by Gerard Ateshian on 11/17/13.
//  Copyright (c) 2013 Columbia University. All rights reserved.
//

#include "FEContinuousFiberDistribution.h"

//-----------------------------------------------------------------------------
FEContinuousFiberDistribution::FEContinuousFiberDistribution(FEModel* pfem) : FEElasticMaterial(pfem)
{
	// set material properties
	AddProperty(&m_pFmat, "fibers"      );
	AddProperty(&m_pFDD , "distribution");
	AddProperty(&m_pFint, "scheme"      );
}

//-----------------------------------------------------------------------------
FEContinuousFiberDistribution::~FEContinuousFiberDistribution() {}

//-----------------------------------------------------------------------------
void FEContinuousFiberDistribution::Init()
{
    // set parent materials
    m_pFmat->SetParent(this);
    m_pFDD->SetParent(this);
    m_pFint->SetParent(this);
    
    // propagate pointers to fiber material and density distribution
    // to fiber integration scheme
    m_pFint->m_pFmat = m_pFmat;
    m_pFint->m_pFDD = m_pFDD;
    
    // initialize fiber integration scheme
    m_pFint->Init();
}
