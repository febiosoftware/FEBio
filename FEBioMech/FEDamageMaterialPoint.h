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

#ifdef WIN32
#define max(a,b) ((a)>(b)?(a):(b))
#endif

class FEDamageCriterion;
class FEDamageCriterionUC;

//-----------------------------------------------------------------------------
// Define a material point that stores the damage variable.
class FEDamageMaterialPoint : public FEMaterialPoint
{
public:
    FEDamageMaterialPoint(FEMaterialPoint *pt) : FEMaterialPoint(pt) { m_pDC = 0; m_pDU = 0; }
    FEDamageMaterialPoint(FEMaterialPoint *pt, FEDamageCriterion* pDC) : FEMaterialPoint(pt) { m_pDC = pDC; m_pDU = 0; }
    FEDamageMaterialPoint(FEMaterialPoint *pt, FEDamageCriterionUC* pDU) : FEMaterialPoint(pt) { m_pDC = 0; m_pDU = pDU; }
    
    FEMaterialPoint* Copy();
    
    void Init(bool bflag);
    
    void ShallowCopy(DumpStream& dmp, bool bsave);
    
    void Serialize(DumpFile& ar);
    
public:
	double	m_Etrial;		//!< trial damage criterion at time t
	double	m_Emax;			//!< max damage criterion up to time t
	double	m_D;			//!< damage (0 = no damage, 1 = complete damage)
    FEDamageCriterion*  m_pDC;      //!< pointer to damage criterion material
    FEDamageCriterionUC*  m_pDU;    //!< pointer to uncoupled damage criterion material
};


#endif /* defined(__FEBioMech__FEDamageMaterialPoint__) */
