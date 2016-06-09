//
//  FEViscousMaterialPoint.hpp
//  FEBioMech
//
//  Created by Gerard Ateshian on 4/7/16.
//  Copyright © 2016 febio.org. All rights reserved.
//

#ifndef FEViscousMaterialPoint_hpp
#define FEViscousMaterialPoint_hpp

#include "FECore/FEMaterialPoint.h"

//-----------------------------------------------------------------------------
// Define a material point that stores the deformation gradient at previous time point.
class FEViscousMaterialPoint : public FEMaterialPoint
{
public:
    FEViscousMaterialPoint(FEMaterialPoint *pt) : FEMaterialPoint(pt) {}
    
    FEMaterialPoint* Copy();
    
    void Init();

    void Update(const FETimeInfo& timeInfo);

    void Serialize(DumpStream& ar);
    
public:
    mat3d VelocityGradient();
    
    mat3ds RateOfDeformation();
    
public:
    mat3d	m_Fp;		//!< deformation gradient at previous time point
};

#endif /* FEViscousMaterialPoint_hpp */
