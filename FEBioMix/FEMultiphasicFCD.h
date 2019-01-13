//
//  FEMultiphasicFCD.hpp
//  FEBioMix
//
//  Created by Gerard Ateshian on 3/8/18.
//  Copyright © 2018 febio.org. All rights reserved.
//
//  This material implements a FEMultiphasicStandard material where
//  inhomogeneous fixed charge density may be specified for each
//  element in the mesh data description.  The FCD at the element
//  level is multiplied by the FCD at the material level, to account
//  for a loadcurve associated with the material level FCD.

#ifndef FEMultiphasicFCD_hpp
#define FEMultiphasicFCD_hpp

#include <FEBioMix/FEMultiphasicStandard.h>

class FECORE_API FEFCDMaterialPoint : public FESolutesMaterialPoint
{
public:
    FEFCDMaterialPoint(FEMaterialPoint* ppt);
    
    void Init(bool bflag);
    
public:
    double    m_cFr;
};

class FEMultiphasicFCD : public FEMultiphasicStandard
{
public:
    FEMultiphasicFCD(FEModel* pfem) : FEMultiphasicStandard(pfem){}
    
    FEMaterialPoint* CreateMaterialPointData() override;
    
    //! fixed charge density
    double FixedChargeDensity(FEMaterialPoint& pt) override;
    
};

#endif /* FEMultiphasicFCD_hpp */
