//
//  FEPrescribedActiveContractionIsotropic.h
//  FEBioMech
//
//  Created by Gerard Ateshian on 3/13/15.
//  Copyright (c) 2015 febio.org. All rights reserved.
//

#ifndef __FEBioMech__FEPrescribedActiveContractionIsotropic__
#define __FEBioMech__FEPrescribedActiveContractionIsotropic__

#include "FEElasticMaterial.h"

//-----------------------------------------------------------------------------
// This material implements an active contraction model which can be used
// as a component of a solid mixture material.
class FEPrescribedActiveContractionIsotropic : public FEElasticMaterial
{
public:
    //! constructor
    FEPrescribedActiveContractionIsotropic(FEModel* pfem);
    
    //! Initialization
    void Init();
    
    //! stress
    mat3ds Stress(FEMaterialPoint& pt);
    
    //! tangent
    tens4ds Tangent(FEMaterialPoint& pt);
    
public:
    double	m_T0;       // prescribed active stress
    
    DECLARE_PARAMETER_LIST();
};


#endif /* defined(__FEBioMech__FEPrescribedActiveContractionIsotropic__) */
