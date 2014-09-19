//
//  FEDamageMaterialPoint.h
//  FEBioMech
//
//  Created by Gerard Ateshian on 9/18/14.
//  Copyright (c) 2014 febio.org. All rights reserved.
//

#ifndef __FEBioMech__FEDamageMaterialPoint__
#define __FEBioMech__FEDamageMaterialPoint__

#include "FECore/FEMaterialPoint.h"

//-----------------------------------------------------------------------------
// Define a material point that stores the damage variable.
class FEDamageMaterialPoint : public FEMaterialPoint
{
public:
	FEDamageMaterialPoint(FEMaterialPoint *pt) : FEMaterialPoint(pt) {}
    
	FEMaterialPoint* Copy()
	{
		FEDamageMaterialPoint* pt = new FEDamageMaterialPoint(*this);
		if (m_pt) pt->m_pt = m_pt->Copy();
		return pt;
	}
    
	void Init(bool bflag)
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
			m_Emax = max(m_Emax, m_Etrial);
		}
        
		// don't forget to intialize the nested data
		if (m_pt) m_pt->Init(bflag);
	}
    
	void ShallowCopy(DumpStream& dmp, bool bsave)
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
    
	void Serialize(DumpFile& ar)
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
    
public:
	double	m_Etrial;		//!< trial damage criterion at time t
	double	m_Emax;			//!< max damage criterion up to time t
	double	m_D;			//!< damage (0 = no damage, 1 = complete damage)
};


#endif /* defined(__FEBioMech__FEDamageMaterialPoint__) */
