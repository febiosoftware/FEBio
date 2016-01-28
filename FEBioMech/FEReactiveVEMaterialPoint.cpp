//
//  FEReactiveVEMaterialPoint.cpp
//  FEBioMech
//
//  Created by Gerard Ateshian on 11/30/14.
//  Copyright (c) 2014 febio.org. All rights reserved.
//

#include "FEReactiveVEMaterialPoint.h"
#include "FEElasticMaterial.h"

///////////////////////////////////////////////////////////////////////////////
//
// FEReactiveVEMaterialPoint
//
///////////////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
//! Create a shallow copy of the material point data
FEMaterialPoint* FEReactiveVEMaterialPoint::Copy()
{
    FEReactiveVEMaterialPoint* pt = new FEReactiveVEMaterialPoint(*this);
    if (m_pNext) pt->m_pNext = m_pNext->Copy();
    return pt;
}

//-----------------------------------------------------------------------------
//! Initializes material point data.
void FEReactiveVEMaterialPoint::Init(bool bflag)
{
    FEElasticMaterialPoint& pt = *m_pNext->ExtractData<FEElasticMaterialPoint>();
    if (bflag)
    {
        // initialize data to zero
        m_Fi.clear();
        m_Ji.clear();
        m_v.clear();
        m_w.clear();
    }
    else
    {
        // check if the current deformation gradient is different from that of
        // the last generation, in which case store the current state
        if (m_pRve) {
            if (m_pRve->NewGeneration(*this)) {
                m_Fi.push_back(pt.m_F.inverse());
                m_Ji.push_back(1./pt.m_J);
                m_v.push_back(FEMaterialPoint::time);
                double w = m_pRve->ReformingBondMassFraction(*this);
                m_w.push_back(w);
            }
        }
        else {
            if (m_pRuc->NewGeneration(*this)) {
                m_Fi.push_back(pt.m_F.inverse());
                m_Ji.push_back(1./pt.m_J);
                m_v.push_back(FEMaterialPoint::time);
                double w = m_pRuc->ReformingBondMassFraction(*this);
                m_w.push_back(w);
            }
        }
    }
    
    // don't forget to initialize the base class
    FEMaterialPoint::Init(bflag);
}

//-----------------------------------------------------------------------------
//! Serialize data to the archive
void FEReactiveVEMaterialPoint::Serialize(DumpStream& ar)
{
    FEMaterialPoint::Serialize(ar);
    
    if (ar.IsSaving())
    {
        int n = (int)m_Fi.size();
        ar << n;
        for (int i=0; i<n; ++i) ar << m_Fi[i] << m_Ji[i] << m_v[i] << m_w[i];
    }
    else
    {
        int n;
        ar >> n;
        for (int i=0; i<n; ++i) ar >> m_Fi[i] >> m_Ji[i] >> m_v[i] >> m_w[i];
    }
}
