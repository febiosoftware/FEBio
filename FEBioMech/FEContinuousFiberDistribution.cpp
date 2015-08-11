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
	m_pFmat.SetName("fibers"      ).SetID(0);
	m_pFDD .SetName("distribution").SetID(1);
	m_pFint.SetName("scheme"      ).SetID(2);
}

//-----------------------------------------------------------------------------
FEContinuousFiberDistribution::~FEContinuousFiberDistribution() {}

//-----------------------------------------------------------------------------
int FEContinuousFiberDistribution::MaterialProperties()
{
	return 3;
}

//-----------------------------------------------------------------------------
//! get a specific material property
FEProperty* FEContinuousFiberDistribution::GetMaterialProperty(int i)
{
	switch(i)
	{
	case 0: return &m_pFmat; break;
	case 1: return &m_pFDD; break;
	case 2: return &m_pFint; break;
	}
	return 0;
}

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
