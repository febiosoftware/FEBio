//
//  FEDamageMaterialPoint.cpp
//  FEBioMech
//
//  Created by Gerard Ateshian on 9/18/14.
//  Copyright (c) 2014 febio.org. All rights reserved.
//

#include "FEDamageMaterialPoint.h"
#include "FEDamageCriterion.h"
#include "FEDamageCriterionUC.h"

FEMaterialPoint* FEDamageMaterialPoint::Copy()
{
    FEDamageMaterialPoint* pt = new FEDamageMaterialPoint(*this);
    if (m_pt) pt->m_pt = m_pt->Copy();
    return pt;
}

void FEDamageMaterialPoint::Init(bool bflag)
{
    if (bflag)
    {
        // intialize data to zero
        m_Emax = 0;
        m_Etrial = 0;
        m_D = 0;
    }
    else
    {
        if (m_pDC) m_Emax = max(m_Emax, m_pDC->DamageCriterion(*this));
        else if (m_pDU) m_Emax = max(m_Emax, m_pDU->DamageCriterion(*this));
        else m_Emax = max(m_Emax, m_Etrial);
    }
    
    // don't forget to intialize the nested data
    if (m_pt) m_pt->Init(bflag);
}

void FEDamageMaterialPoint::ShallowCopy(DumpStream& dmp, bool bsave)
{
    if (bsave)
    {
        dmp << m_Etrial << m_Emax << m_D;
    }
    else
    {
        dmp >> m_Etrial >> m_Emax >> m_D;
    }
    if (m_pt) m_pt->ShallowCopy(dmp, bsave);
}

void FEDamageMaterialPoint::Serialize(DumpFile& ar)
{
    if (ar.IsSaving())
    {
        ar << m_Emax;
    }
    else
    {
        ar >> m_Emax;
    }
}
