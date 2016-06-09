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

//-----------------------------------------------------------------------------
// Define a material point that stores the damage variable.
class FEDamageMaterialPoint : public FEMaterialPoint
{
public:
    FEDamageMaterialPoint(FEMaterialPoint *pt) : FEMaterialPoint(pt) {}
    
    FEMaterialPoint* Copy();
    
    void Init();
    void Update(const FETimeInfo& timeInfo);
    
    void Serialize(DumpStream& ar);
    
public:
	double	m_Etrial;		//!< trial damage criterion at time t
	double	m_Emax;			//!< max damage criterion up to time t
	double	m_D;			//!< damage (0 = no damage, 1 = complete damage)
};


#endif /* defined(__FEBioMech__FEDamageMaterialPoint__) */
