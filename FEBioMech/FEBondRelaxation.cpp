//
//  FEBondRelaxation.cpp
//  FEBioMech
//
//  Created by Gerard Ateshian on 8/25/14.
//  Copyright (c) 2014 febio.org. All rights reserved.
//

#include "FEBondRelaxation.h"

//-----------------------------------------------------------------------------
// Material parameters for FEBondRelaxation
void FEBondRelaxation::Init()
{
	FEMaterial::Init();
}

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_PARAMETER_LIST(FEBondRelaxationConstIso, FEBondRelaxation)
    ADD_PARAMETER(m_t[0], FE_PARAM_DOUBLE, "t1");
    ADD_PARAMETER(m_t[1], FE_PARAM_DOUBLE, "t2");
    ADD_PARAMETER(m_t[2], FE_PARAM_DOUBLE, "t3");
    ADD_PARAMETER(m_t[3], FE_PARAM_DOUBLE, "t4");
    ADD_PARAMETER(m_t[4], FE_PARAM_DOUBLE, "t5");
    ADD_PARAMETER(m_t[5], FE_PARAM_DOUBLE, "t6");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! Constructor.
FEBondRelaxationConstIso::FEBondRelaxationConstIso(FEModel* pfem) : FEBondRelaxation(pfem)
{
    m_nt = 0;
	for (int i=0; i<MAX_TERMS; ++i) m_t[i] = 0;
}

//-----------------------------------------------------------------------------
//! Initialization.
void FEBondRelaxationConstIso::Init()
{
    FEBondRelaxation::Init();
    
    m_nt = 0;
	for (int i=0; i<MAX_TERMS; ++i) {
        if (m_t[i] < 0) throw MaterialError("relaxation times must be > 0");
        if (m_t[i] > 0) ++m_nt;
    }
    if (m_nt == 0)
        throw MaterialError("at least one relaxation time must be > 0");
}

//-----------------------------------------------------------------------------
//! Relaxation function
double FEBondRelaxationConstIso::Relaxation(FEMaterialPoint& mp, const double t)
{
	// --- constant relaxation times ---
    double g = 0;
    
    for (int i=0; i<m_nt; ++i)
        if (m_t[i] > 0) {
            g += exp(-t/m_t[i]);
        }
	
	return g/m_nt;
}
