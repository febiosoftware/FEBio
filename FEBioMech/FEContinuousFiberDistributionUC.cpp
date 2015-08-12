//
//  FEContinuousFiberDistributionUC.cpp
//  FEBioMech
//
//  Created by Gerard Ateshian on 8/5/14.
//  Copyright (c) 2014 febio.org. All rights reserved.
//

#include "FEContinuousFiberDistributionUC.h"

//-----------------------------------------------------------------------------
FEContinuousFiberDistributionUC::FEContinuousFiberDistributionUC(FEModel* pfem) : FEUncoupledMaterial(pfem)
{
	// set material properties
	AddProperty(&m_pFmat, "fibers"      );
	AddProperty(&m_pFDD , "distribution");
	AddProperty(&m_pFint, "scheme"      );
}

//-----------------------------------------------------------------------------
FEContinuousFiberDistributionUC::~FEContinuousFiberDistributionUC() {}

//-----------------------------------------------------------------------------
void FEContinuousFiberDistributionUC::Init()
{
    m_K = m_pFmat->m_K;
    
	// propagate pointers to fiber material and density distribution
    // to fiber integration scheme
    m_pFint->m_pFmat = m_pFmat;
    m_pFint->m_pFDD = m_pFDD;
    
    // initialize fiber integration scheme
    FEUncoupledMaterial::Init();
}
